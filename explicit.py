import numpy as np
import math
#Liner properties
k_insulator = 0.2                                          #Insulation thermal conductivity [W/m-K]
rho_insulator = 860                                        #Insulation density [kg/m3]
cp_insulator = 2000                                        #Insulation specific heat capacity [J/kg-K]
alpha_insulator = float(k_insulator / (cp_insulator * rho_insulator))  #Insulation thermal diffusivity [m^2/s]
#------------------------------------------------
#Case properties
k_case = 167                                      #Casing thermal conductivity [W/m-K]
rho_case = 2700                                   #Casing density [kg/m3]
cp_case = 896                                     #Casing specific heat capacity [J/kg-K]
alpha_case = float(k_case / (cp_case * rho_case)) #Casing thermal diffusivity [m^2/s]

#-----------------------------------------------
#Convection Coefficient
h_m = 1295  #Average convection coefficient [W/m^2-K]

#------------------------------------------------
#Inputs
#Radial coordinate
insulator_thk = 0.003         #Insulation thickness [m]
case_thk = 0.00538            #Case thickness |VERIFY THIS VALUEEEE| [m]
r_1 = 0.11                    #Radial position of insulation beginning [m]
r_2 = r_1 + insulator_thk     #Radial position of insulation/casing interface [m]
r_3 = r_2 + case_thk          #Radial position of casing end [m]
sectionsr = 100                #Nº of sections between r_1 and r_3
dr = (r_3 - r_1) / sectionsr  #r variation [m]
nr = sectionsr + 1            #Nº of radial points

#Height coordinate
z = 1.42         #Height of combustion chamber [m]
sectionsz = 25   #Nº of sections beewtween 0 and z
dz = z/sectionsz #z variation [m]
nz = sectionsz+1 #Nº of z points

#Time
t = 7                                                                        #Burn time [s]
dt_insulator = (k_insulator*dr**2)/(4*alpha_insulator*(k_insulator+h_m*dr))  #Stability dt for insulator [s]
dt_case = (k_case*dr**2)/(4*alpha_case*(k_case+h_m*dr))                      #Stability dt for casing [s]
if dt_insulator >= dt_case: dt = dt_case
else: dt = dt_insulator                                                      #t variation [s]

nt = math.ceil(t / dt)                                                       #Nº of time points

print(nt)

d_insulator = alpha_insulator * (dt / (dr**2))  #Insulation Fourier number
d_case = alpha_case * (dt / (dr**2))            #Casing Fourier number

#Temperature
Ta = 297   #Ambient temperature in Tatuí-SP-Brazil in the month of august[K]
Tc = 1600  #Chamber temperature [K] -> This is a given value, yet this is not the real value

#we have to set the initial/boundary contitions
r = [r_1 + i * dr for i in range(nr)]

T = np.zeros([nt, nr, nz])
for i in range(nr):
  for j in range(nz):
    T[0][i][j] = Ta

#Calculating the T-matrix
for j in range(nt - 1):
    for i in range(nr):
      for l in range(nz):
        if r[i] <= r_2:  #Set insulation material
            alpha = alpha_insulator
            rho = rho_insulator
            cp = cp_insulator
            k = k_insulator
        if r[i] > r_2:   #Set casing material
            alpha = alpha_case
            rho = rho_case
            cp = cp_case
            k = k_case
        if i == nr - 1:
          if l == nz - 1:   #r end and z end
            T[j + 1][i][l] = T[j][i][l] + alpha * dt * (((T[j][i][l]-2*T[j][i][l-1]+T[j][i][l-2])/(dz**2))+((T[j][i][l] - T[j][i - 1][l])/(r[i] * dr)) + ((T[j][i][l] + T[j][i-2][l] - 2 * T[j][i-1][l]) / (dr**2)))
          elif l == 0:      #r end and z start
            T[j + 1][i][l] = T[j][i][l] + alpha * dt * (((T[j][i][l]-2*T[j][i][l+1]+T[j][i][l+2])/(dz**2))+((T[j][i][l] - T[j][i - 1][l])/(r[i] * dr)) + ((T[j][i][l] + T[j][i-2][l] - 2 * T[j][i-1][l]) / (dr**2)))
          else:             #r end and z middle
            T[j + 1][i][l] = T[j][i][l] + alpha * dt * (((T[j][i][l+1]-2*T[j][i][l]+T[j][i][l-1])/(dz**2))+((T[j][i][l] - T[j][i - 1][l])/(r[i] * dr)) + ((T[j][i][l] + T[j][i-2][l] - 2 * T[j][i-1][l]) / (dr**2)))
        elif i == 0:
          if l == nz-1:     #r start and z end
            T[j + 1][i][l] = T[j][i][l] + ((2 * dt) / (rho * cp * dr * dz)) * (h_m * dz * (Tc - T[j][i][l]) + k * dz * (T[j][i+1][l] - T[j][i][l])/dr + k * dr * (T[j][i][l] -T[j][i][l-1])/(2 * dz))
          elif l == 0:      #r start and z start
            T[j + 1][i][l] = T[j][i][l] + ((2 * dt) / (rho * cp * dr * dz)) * (h_m * dz * (Tc - T[j][i][l]) + k * dz * (T[j][i+1][l] - T[j][i][l])/dr + k * dr * (T[j][i][l+1] - T[j][i][l])/(2 * dz))
          else:             #r start and z middle
            T[j + 1][i][l] = T[j][i][l] + ((2 * dt) / (rho * cp * dr * dz)) * (h_m * dz * (Tc - T[j][i][l]) + k * dz * (T[j][i+1][l] - T[j][i][l])/dr + k * dr * (T[j][i][l+1] - T[j][i][l-1])/(2 * dz))
        else:
          if l == nz-1:     #r middle and z end
            T[j + 1][i][l] = T[j][i][l] + alpha * dt * ((
                (T[j][i+1][l] - T[j][i-1][l])/(2 * r[i] * dr)) + ((T[j][i+1][l] + T[j][i-1][l-2] - 2 * T[j][i][l]) / (dr**2)) + ((T[j][i][l] + T[j][i][l-2] - 2*T[j][i][l-1]) / (dz**2)))
          elif l == 0:      #r middle and z start
            T[j + 1][i][l] = T[j][i][l] + alpha * dt * ((
                (T[j][i+1][l] - T[j][i-1][l])/(2 * r[i] * dr)) + ((T[j][i+1][l] + T[j][i-1][l-2] - 2 * T[j][i][l]) / (dr**2)) + ((T[j][i][l] + T[j][i][l+2] - 2*T[j][i][l+1]) / (dz**2)))
          else:             #r middle and z middle
            T[j + 1][i][l] = T[j][i][l] + alpha * dt * ((
                (T[j][i+1][l] - T[j][i-1][l])/(2 * r[i] * dr)) + ((T[j][i+1][l] + T[j][i-1][l-2] - 2 * T[j][i][l]) / (dr**2)) + ((T[j][i][l+1] + T[j][i][l-1] - 2*T[j][i][l]) / (dz**2)))

print(T)
