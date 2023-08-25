import numpy as np

# constant
Q = 1.6e-19  # [C]電荷（クーロン単位）
EPS = 8.9e-12  # [F m*-1]真空の誘電率
MYU = 1.3e-6  # [N A**-2]真空の有事率
ME = 9.109e-31  # [kg]電子の質量
MH = 1.7e-27  # [kg]水素イオンの質量
MHE = 6.65e-27  # [kg]ヘリウムイオンの質量
MO = 2.7e-26  # [kg]酸素イオンの質量
C = 2.97e8  # [m s**-1]光速
L = 6.0  # L-value
ML = 0 # 磁気緯度
B_E = 3.11*10**-5  # [T]地球の表面、磁気赤道上の磁場の大きさ
B0 = B_E/L**3*(1-3*np.sin(ML)**2)**(1/2)/np.cos(ML)**6  # [T]磁場の大きさ

# plasma parameter
NE = 2e6  # [m**-3]電子密度
ion_ratio = np.array([0.24, 0.08, 0.68]) # 各イオン種の電子密度に対する割合
nh = ion_ratio[0]*NE # 水素イオンの密度
nhe = ion_ratio[1]*NE # ヘリウムイオンの密度
no = ion_ratio[2]*NE # 酸素分子イオンの密度

RHO = sum(ion_ratio*np.array([MH, MHE, MO]))
# 電子密度に対するイオン種の割合にそれぞれのイオン質量をかけて総和をとる

pi_e = {'value': (NE*Q**2/(EPS*ME))**0.5, 'label': r'$\Pi_{e}$'}
pi_h = {'value': (nh*Q**2/(EPS*MH))**0.5, 'lable': r'$\Pi_{h}$'}
pi_he = {'value': (nhe*Q**2/(EPS*MHE))**0.5, 'label': r'$\Pi_{he}$'}
pi_o = {'value': (no*Q**2/(EPS*MO))**0.5, 'label': r'$\Pi_{o}$'}
# プラズマ周波数

omega_e = {'value': -Q*B0/ME, 'label': r'$\Omega_{ce}$'}
omega_h = {'value': Q*B0/MH, 'label': r'$\Omega_{ch}$'}
omega_o = {'value': Q*B0/MO, 'label': r'$\Omega_{co}$'}
omega_he = {'value': Q*B0/MHE, 'label': r'$\Omega_{che}$'}
# サイクロトロン周波数

wlh = {'value': np.sqrt((omega_h['value']**2 + pi_h['value']**2) / (1 + (pi_e['value']/omega_e['value'])**2)),
       'label': r'$\Omega_{lh}$'}
wuh = {'value': np.sqrt(omega_e['value']**2 + pi_e['value']**2),
       'label': r'$\Omega_{uh}$'}
# 定義に関しては4章を参照。垂直方向に伝搬する電磁波の共鳴周波数。
# UH周波数は、電子のサイクロトロン周波数よりも大きく、イオンのサイクロトロン周波数よりはるかに大きいので電子のみで計算している。
# LH周波数はイオンのサイクロトロン周波数と電子のプラズマ周波数の間にある。ここでは水素イオンと電子のみで計算している。

Va = B0/(MYU*RHO)**0.5