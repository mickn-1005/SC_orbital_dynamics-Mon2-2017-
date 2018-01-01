# -*- coding: utf-8 -*-
"""
Created on Tue May 30 08:51:16 2017

システム学第２最終レポート
"""

#pythonでNORADのTLEを読み込んで軌道計算。

import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.pyplot as mpl
import importlib
importlib.import_module('mpl_toolkits.mplot3d').__path__

from mpl_toolkits.mplot3d import Axes3D


def split_tle(code):
    """
    ex.環境観測衛星「みどり」

    MIDORI(ADEOS)
    1 24277U 96046A   09116.47337938 -.00000023  00000-0  73445-5 0   432
    2 24277  98.3597  83.2073 0002090  64.7512 295.3886 14.28595439661547
    """
    #elementへの分解
    print("\nshaping elements\n")
    elements = code.split()
    for i in range(len(elements)):
        print(elements[i])  # 分割できていることがわかった。衛星名をエレメントから消すか軌道情報だけを半βつする仕組みが必要。
    print("\npop exceed infos...")
    for i in range(len(elements)):
        if(elements[i] == "1"):
            if(("S" in elements[i+1])or("U" in elements[i+1])):  #TLEの軌道情報の最初
                if(elements[i+9] in elements[i+1]):  #絞り込む
                    for j in range(i):
                        elements.pop(0)
                    print("info done!!")
                    break
        i += 1
    # print("new list is\n")
    # print(elements)  #分割できている!
    return elements

    """
    elementsの各軌道要素
    00;行数（１）
    01;衛星通し番号と公開か否か
    02;国際衛星識別符号
    03;軌道元期（２桁年＋日単位。最小桁の制度は0.864ms）
    04;平均運動の１時微分値の２分の１（回転/day^2）
    05;２次微分の６分の１（回転/day^3）使用されない
    06;B*抗力項下２桁は１０の指数
    07;軌道モデル種別
    08;軌道要素通番+チェックサム（下１桁）
    09;行数（２）
    10;カタログ番号（01から"S"または"U"を消去したもの）
    11;軌道傾斜角
    12;昇交点赤経
    13;離心率
    14;近地点引数
    15;平均近点角
    16;平均運動＋通算周期数（４桁）＋チェックサム（下１桁）
    """

def get_height(time, line):
    """
    軌道面内の衛星高度を取得
    """
    GM = 2.975537 * pow(10,15)  #地心重力定数
    Mm = float(line[16][:12]) + float(line[4])*2 * time #time経過後の１次運動
    # print(Mm)
    kep3 = GM / (4*pow(math.pi,2)*pow(Mm,2))
    r_long = pow(kep3, 1/3)  #軌道長半径を取得できた
    height = r_long - r_earth  #赤道半径を引く やや精度不足？？
    return height, r_long #高度を取得できた

def get_eccen_anom(time, line):
    """
    ニュートン法によって離心近点角を取得
    """
    M = float(line[15]) / 360 + float(line[16][:12]) * time + float(line[4])*0.1*pow(time,2) # 平均近点角を取得した単位（rev)\
    # print("average anom is")
    # print(M)
    e = float(line[13])*0.00000001
    sat_angle, sat_dev = math.modf(M)
    sat_angle = sat_angle * 360  #平均近点角を取得空いた単位（°）
    #離心率とMから離心近点角を取得。初期値M
    def ecce_f(x):
        value = x - e * np.sin(x) - sat_angle
        return value

    def new_rap(x):
        value = x - (x - e * np.sin(x) - sat_angle) / (1 - e * np.cos(x))
        return value

    M_ec = sat_angle #初期値
    # print(ecce_f(M_ec))
    while abs(ecce_f(M_ec))>0.0000001:
        # print("ecce_f is")
        # print(ecce_f(M_ec))
        M_ec = new_rap(M_ec)
    # print('eccemtric anomaly is')
    # print(M_ec)  # 離心近点角を取得できた

    return M_ec

# 3次元回転行列の取得
def rot_3Dmatz(theta):
    co = np.cos(theta)
    si = np.sin(theta)

    matrix_3 = np.array([[co, -si, 0], [si, co, 0], [0,0,1]])

    return matrix_3

def rot_3Dmatx(theta):
    co = np.cos(theta)
    si = np.sin(theta)

    matrix_3 = np.array([[1,0,0],[0,co,-si],[0,si,co]])

    return matrix_3


def get_earth_position( time,line, eccen_anom, maj_axis ):
    """
    地心３次元座標内の衛星位置の取得
    """
    #軌道面内の衛星の座標を取得
    omega0 = float(line[14])  #近地点引数
    i = float(line[11])  #軌道傾斜角
    omega_cap0 = float(line[12])  #昇交点赤経
    e = float(line[13])*0.00000001 #離心率

    U = maj_axis * np.cos(eccen_anom) - e * maj_axis
    V = maj_axis * np.sqrt(1 - pow(e, 2)) * np.sin(eccen_anom)
    orb_pos = np.array([[U],[V],[0]])#軌道面内座標取得
    # 軌道要素の補正
    omega = omega0 + time * 180 * 0.174 * (2 - 2.5*pow(np.sin(i),2)) / (math.pi * pow(maj_axis / r_earth, 3.5))
    omega_cap = omega_cap0 - time * 180 * 0.174 * np.cos(i) / (math.pi * pow(maj_axis / r_earth, 3.5))

    #地心３次元座標への変換
    total_rot = np.dot(np.dot(rot_3Dmatz(omega_cap),rot_3Dmatx(i)), rot_3Dmatz(omega))
    ear_pos = np.dot(total_rot, orb_pos)
    # print(ear_pos)  #地心座標の取得
    return ear_pos

def plot_sphere():
    """
    地球表面の生成
    """
    # Make data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r_earth * np.outer(np.cos(u), np.sin(v))
    y = r_earth * np.outer(np.sin(u), np.sin(v))
    z = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='b')




#メイン関数
tle = "MIDORI (ADEOS)\n1 24277U 96046A   09116.47337938 -.00000023  00000-0  73445-5 0   432\n2 24277  98.3597  83.2073 0002090  64.7512 295.3886 14.28595439661547"
mod_tle = split_tle(tle)

r_earth = 6378.137 # 地球の赤道半径単位：km
year = 2017
month = 6
day = 15  #とりあえず日単位で出せるように
r_earth = 6378.1366 # 地球の赤道半径

#matplotlib関係の宣言
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_zlabel('NS axis')

delta = 7.18  #初期値（元期起点）単位：day
delta_t = 0.000001 #時間差単位：day
m_cnt = 0 #試行回数
b_count = 100000 #break値, 42442においてノイズ発生
x = np.zeros(b_count+2)
y = np.zeros(b_count+2)
z = np.zeros(b_count+2)

"""以下、mainの処理"""
print("start calculation")
while True:
    sat_h, arm = get_height(delta, mod_tle)
    anom = get_eccen_anom(delta, mod_tle)
    position = get_earth_position(delta, mod_tle,anom,arm)
    x[m_cnt] = position[0]
    y[m_cnt] = position[1]
    z[m_cnt] = position[2]
    if m_cnt > b_count:
        print("break from m_cnt")
        break
    delta += delta_t
    m_cnt += 1

print("calc done.\nplotting...")

"""plotting"""
ax.scatter3D(x[0],y[0],z[0],'^', color = 'black' )
ax.plot(x,y,z, color = 'red')
plot_sphere()
plt.show()
