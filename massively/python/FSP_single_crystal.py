import numpy as np
pi = np.pi

### NMR parameters ############################################################
I          = 3/2 #核スピン（分数でも少数でもOK）
gyr        = 11.289 # 核スピン磁気回転比（MHz/T）
H_0        = 10.0 # 磁場（T）
K_shift    = 0.000 / 100 # ナイトシフト（%）
nu_Q       = 9.28 # NQR共鳴周波（MHz）
eta        = 0.000 # 非対称パラメータ
theta      = 20.0 * pi / 180.0 # 電場勾配の最大主軸をZ方向とする.
phi        = 00.0 * pi / 180.0 # 磁場方向を極座標(theta,phi)で表す．単位はrad.
alpha      = 00.0 * pi / 180.0 # z軸をx軸に向ける回転方向（オイラー角. 単位はrad.）
beta       = 00.0 * pi / 180.0 # y軸をz軸に向ける回転方向（結晶主軸(x,y,z)と電場勾配の主軸(X,Y,Z)のズレを補正する.） 
gamma      = 00.0 * pi / 180.0 # x軸をy軸に向ける回転方向
N          = 1024 # 刻み幅
sigma      = 0.2 # 標準偏差（MHz）.Intensityは"W[i,j] * np.exp(-(x - d_f[i,j])**2 / (2*sigma**2))"の重ね合わせ. xの刻みより大きくする.
### End of parameters #########################################################

nu_Q = nu_Q / (1 + (eta**(2)/3))**(1/2)
# I=3/2でetaが有限のとき，ゼロ磁場NQRにおける1本のスペクトルの共鳴周波数を固定するためにnu_Qを再設定.

### Nuclear spin ##############################################################
_2Ip1    = int(2 * I + 1)
I_2     = I * (I + 1)
m      = np.arange(I, -I - 1.0, -1.0)

I_plus = np.diag(np.sqrt(I_2 - (m[0:-1]) * (m[0:-1] - 1.0)), 1)
I_minus = np.conj(I_plus.T)
I_x = +0.5  * (I_plus + I_minus)
I_y = -0.5j * (I_plus - I_minus)
I_z = np.diag(m)

### Rotaion matric ############################################################
R_x_beta  = np.array([[1, 0, 0], [0, np.cos(beta), -np.sin(beta)],[0,np.sin(beta),np.cos(beta)]])
R_y_alpha = np.array([[np.cos(alpha), 0, np.sin(alpha)], [0, 1, 0],[-np.sin(alpha), 0, np.cos(alpha)]])
R_z_gamma = np.array([[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0],[0, 0, 1]])
R_xy    = np.dot(R_x_beta, R_y_alpha)
R_zxy   = np.dot(R_z_gamma, R_xy)

### Coordinate transformation #################################################
I_X     = I_x * R_zxy[0,0] + I_y * R_zxy[0,1] + I_z * R_zxy[0,2] 
I_Y     = I_x * R_zxy[1,0] + I_y * R_zxy[1,1] + I_z * R_zxy[1,2] 
I_Z     = I_x * R_zxy[2,0] + I_y * R_zxy[2,1] + I_z * R_zxy[2,2] 
I_Plus  = I_X + 1.0j * I_Y
I_Minus = I_X - 1.0j * I_Y

### Hamiltonian ###############################################################
Hamiltonian = -gyr * H_0 * (1 + K_shift) * (
        I_z * np.cos(theta) + (I_x * np.cos(phi) + I_y * np.sin(phi)) * np.sin(theta)
        ) + nu_Q / 6.0 * (3 * np.dot(I_Z, I_Z) - I_2 * np.eye(_2Ip1)
                         + 0.5 * eta * (np.dot(I_Plus, I_Plus)+np.dot(I_Minus, I_Minus)))

# hbarは省力しているので，固有値の次元は周波数[MHz]となる.
# ここでH_0とIの内積を"I_z * np.cos(theta) + (I_x * np.cos(phi) + I_y * np.sin(phi)) * np.sin(theta)"としている.
# すなわち結晶主軸(x,y,z)基準で磁場H_0の方向を極座標(theta, phi)で表している.
# 電場勾配の最大主軸をZ方向としているので, 電気四重極ハミルトニアンはI_Z, I_Plus, I_Minusを用いている.

### Eigenvector and energy eigenvalue ##########################################
f, v = np.linalg.eigh(Hamiltonian, 'L')

zero  = np.zeros((_2Ip1, _2Ip1, ))
f_m   = np.add(f, zero)
d_f    = np.abs(f_m.T - f_m)

### Transition probability #####################################################
w_plus  = np.dot(np.dot( np.conj(v.T),I_Plus),v)
w_minus = np.dot(np.dot( np.conj(v.T),I_Minus),v)
W = abs(w_plus * w_plus) + abs(w_minus * w_minus)

###Gaussian####################################################################
x       = np.arange(0, 2*np.amax(d_f), np.amax(d_f)/N) # 0～2*np.amax(d_f)までnp.amax(d_f)/N刻みの数値の配列.

import matplotlib.pyplot as plt

Intensity = 0
for i in range(int(2*I+1)):
    for j in range(int(2*I+1)):
        if i==j:
            pass
        else:
            Intensity = Intensity + W[i,j] * np.exp(-(x - d_f[i,j])**2 / (2*sigma**2))

Intensity = Intensity / np.amax(Intensity) # np.amax(Intensity)で規格化.
# a = 1.2 #強度の倍率
# Intensity = a*Intensity

print("Energy eigenvalue")
print(f)
print("Eigenwavefunction")
print(v)
print("Resonance frequency")
print(d_f)
print("Transition probability")
print(W)

fig = plt.figure(figsize = (8, 6))
ax = fig.add_subplot(111)

ax.set_xlabel("freq [MHz]", fontsize = 14)
ax.set_ylabel("Intensity", fontsize = 14)

d_f_non0 = d_f *W /(W+0.01) # W[i,j]=0のd_f[i,j]を消去した行列を用意.
ax.set_xlim([0, np.amax(d_f_non0)*1.2]) # 軸範囲を設定.
ax.set_ylim([0, 1.2])
print(d_f_non0)

ax.plot(x, Intensity, color = "red", label = "Intensity")
ax.legend(fontsize = 14)
plt.show()