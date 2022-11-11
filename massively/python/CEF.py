import numpy as np

###Calculation Stevens Operators#################################################

###Select Unit###################################################################
Unit = 1 # 1:emu/mol, 2:μ_B/T  , emu/g 
M    = 1833.14 #Molecular weight [g]
R    = 8.31447 #J・K^-1 mol^-1

###J,L,S#########################################################################
J = 4 #整数または分数で入力する.
L = 5
S = 1

n_4f = 1

N_A  = 6.022136 * 10**23 #1/mol
mu_B = 9.2740154 * 10**(-24) #J/T
k_B  = 1.380658 * 10**(-23) #J/K

g_J = 1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))

f4_int = [0,0,12,15,60,14,60,60,60] #Jが整数の場合のF(4)
f4_half = [0,0,0,0,0,60,0,60,0,60,0,60,0,60,0,60]#Jが半整数の場合のF(4)

f6_int = [0,0,0,180,1260,1260,7560,3780]#Jが整数の場合のF(6)
f6_half = [0,0,0,0,0,0,0,1260,0,2520,0,3780,0,7560,0,13860] #Jが半整数の場合のF(6)

if isinstance(J, int):
    f4 = f4_int[J]
    f6 = f6_int[J]
else:
    f4 = f4_half[round(J*2)]
    f6 = f6_half[round(J*2)]

print(g_J)

###Matrix#########################################################################
_2Jp1  = int(2 * J + 1)
J_2    = J * (J + 1)
m      = np.arange(J, -J - 1.0, -1.0)

E1 = np.eye(_2Jp1)
J_z =  np.diag(m)
J_plus = np.diag(np.sqrt(J_2 - (m[0:-1]) * (m[0:-1] - 1.0)), 1)
J_minus = np.conj(J_plus.T)
J_x = +0.5  * (J_plus + J_minus)
J_y = -0.5j * (J_plus - J_minus)
J_2_mat = J_2*E1

#print(J_z)
#print(J_plus)
#print(J_minus)
#print(J_x)
#print(J_y)

O_20 = 3*J_z@J_z - J_2_mat
O_22 = (J_plus@J_minus + J_minus@J_plus)/2
O_xy = (J_x@J_y + J_y@J_x)/2
O_yz = (J_z@J_x + J_x@J_y)/2
O_zx = (J_z@J_x + J_x@J_z)/2
O_40 = 35*J_z@J_z@J_z@J_z - 30*J_2_mat@J_z@J_z + 25*J_z@J_z - 6*J_2_mat + 3*J_2_mat@J_2_mat
O_42 = ( (J_plus@J_plus + J_minus@J_minus) @ (7*J_z@J_z - J_2_mat -5*E1 ) + (7*J_z@J_z - J_2_mat - 5*E1) @ (J_plus@J_plus + J_minus@J_minus) )/4
O_43 = ( J_z@(J_plus@J_plus@J_plus + J_minus@J_minus@J_minus) + (J_plus@J_plus@J_plus + J_minus@J_minus@J_minus)@J_z )/4
O_44 = (J_plus@J_plus@J_plus@J_plus + J_minus@J_minus@J_minus@J_minus)/2
O_60 = 231*J_z@J_z@J_z@J_z@J_z@J_z - (315*J_2_mat - 753*E1) @J_z@J_z@J_z@J_z + (105*J_2_mat@J_2_mat - 525*J_2_mat+294*E1)@J_z@J_z -5*J_2_mat@J_2_mat@J_2_mat + 40*J_2_mat@J_2_mat - 60*J_2_mat

O_62 = ( (33*J_z@J_z@J_z@J_z - 18* J_z@J_z@J_2_mat - 123*J_z@J_z + 10*J_2_mat + 102*E1) @ (J_plus@J_plus + J_minus@J_minus) + (J_plus@J_plus + J_minus@J_minus) @ ( 33*J_z@J_z@J_z@J_z - 18* J_z@J_z@J_2_mat - 123*J_z@J_z + 10*J_2_mat + 102*E1))/4
O_63 = ( (11*J_z@J_z@J_z - 3* J_z@J_2_mat - 59*J_z)@(J_plus@J_plus@J_plus + J_minus@J_minus@J_minus) + (J_plus@J_plus@J_plus + J_minus@J_minus@J_minus) @ (11*J_z@J_z@J_z - 3* J_z@J_2_mat - 59*J_z) )/4
O_64 = ( (J_plus@J_plus@J_plus@J_plus + J_minus@J_minus@J_minus@J_minus)@(11*J_z@J_z - J_2_mat -38*E1) + (11*J_z@J_z - J_2_mat -38*E1)@(J_plus@J_plus@J_plus@J_plus + J_minus@J_minus@J_minus@J_minus) )/4
O_66 = ((J_plus@J_plus@J_plus@J_plus@J_plus@J_plus - J_minus@J_minus@J_minus@J_minus@J_minus@J_minus)/2)

# cubic(立方), hexagonal(六方), tetragonal(正方)ではO_43, O_63が出てこない. 
# O_60はSympyやMathematicaの計算結果と合わない...間違っているかも? 精度の問題?

#print(O_40)
#print(O_44)
#print(O_60)
#print(O_64) 

### Calculation of CEF Hamiltonian #############################################

#(selection CEF Hamiltonian type)  1:Standard Style, 2: LLW Style
S = 2

# Input CEF Parameters
B_20 = 0.01
B_22 = 0.0
B_40 = 0.01
B_42 = 0.0
B_44 = 0.01
B_60 = 0.01
B_62 = 0.0
B_64 = 0.0
B_66 = 0.01

# Input CEF Parameters
WW = -5.3
XX = -0.8

if S == 1:
    H_CEF = B_20*O_20 + B_22*O_22 + B_40*O_40 + B_42*O_42 + B_44*O_44 + B_60*O_60 + B_62*O_62 + B_64*O_64 + B_66*O_66
elif f6 == 0:
    H_CEF = WW*( (XX/f4) *(O_40+5*O_44) )
else:
    H_CEF = WW*( (XX/f4) *(O_40+5*O_44) + (1-abs(XX))*(O_60 - 21*O_64)/f6 )

print(H_CEF)

E_CEF, psi_CEF = np.linalg.eigh(H_CEF, 'L')
E_CEF =  E_CEF - min(E_CEF)

#print(psi_CEF)
print(E_CEF)

J_z_CEF = np.dot(np.dot( np.conj(psi_CEF.T),J_z),psi_CEF)

#print(J_z_CEF)