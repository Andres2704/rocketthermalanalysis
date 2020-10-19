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
radius_outer = 0.07           # Chamber outer radius [m]
radius_inner = 0.064          # Chamber inner radius [m]
sectionsr = 100              #Nº of sections between r_1 and r_3
nr = sectionsr + 1
dr = (radius_outer)/ sectionsr  #r variation [m]
r_1 = dr                 #Radial position of insulation beginning [m]
r_2 = r_1 + radius_inner      #Radial position of insulation/casing interface [m]
r_3 = radius_outer          #Radial position of casing end [m]


#Height coordinate
insulator_thk = 0.003         #Insulation thickness [m]
z = 0.007        #Thickness of bulkhead [m]
z_1 = 0
z_2 = insulator_thk
z_3 = z_2 + z
dz = dr/1.4 #z variation [m]
nz = math.ceil(z/dz) #Nº of z points

#Time
t = 7                                                                  #Burn time [s]

dt = ((dr**2)*(dz**2))/(1.99*alpha_case*(dr**2+dz**2))

nt = math.ceil(t / dt)                                                       #Nº of time points

d_insulator = alpha_insulator * (dt / (dr**2))  #Insulation Fourier number
d_case = alpha_case * (dt / (dr**2))            #Casing Fourier number

#Temperature
Ta = 297   #Ambient temperature in Tatuí-SP-Brazil in the month of august[K]
Tc = 1600  #Chamber temperature [K] -> This is a given value, yet this is not the real value

#we have to set the initial/boundary contitions
r = [r_1 + i * dr for i in range(nr)]
z = [z_1 + i * dr for i in range(nz)]

"""""
print(dr)
print(dz)
print(nz)

print("Insulator")
print("Middle z and start r:")
print((dz**2)/(2*alpha_case))
print("Middle z and middle r:")
print(1/((2*alpha_case/(dr**2)+2*alpha_case/(dz**2))))
print(((dr**2)*(dz**2))/(2*alpha_case*(dr**2+dz**2)))
print("Middle z and end r:")
print((dz**2)*(dr**2)/(2*alpha_case*(dr**2-dz**2)))
print("Start z and start r")
print((rho_case*cp_case)/(k_case/dr**2+k_case/dz**2+2*h_m/dz))
print("Start z and middle r")
print((rho_case*cp_case*dz**2)/(k_case+2*h_m*dz))
print("End z and Middle r")
print((dr**2)*(dz**2)/(alpha_case*(2*dz**2-dr**2)))
"""""

T = np.zeros([nt, nr, nz])
for i in range(nr):
  for j in range(nz):
    T[0][i][j] = Ta

#Calculating the T-matrix
for i in range(nt - 1):
    for l in range(nz):
      for j in range(nr):
        if z[l] <= z_2:  #Set insulation material
            alpha = alpha_insulator
            rho = rho_insulator
            cp = cp_insulator
            k = k_insulator
        if z[l] > z_2:   #Set casing material
            alpha = alpha_case
            rho = rho_case
            cp = cp_case
            k = k_case
        if l == nz-1:
            if j == nr-1:   #z end and r end
                T[i + 1][j][l] = (1 + (alpha * dt) / (r[j] * dr) + (alpha * dt) / (dr ** 2) + (alpha * dt) / (dz ** 2)) * T[i][j][l] + ((-2*alpha*dt)/(dr**2) - (alpha*dt)/(r[j]*dr))*T[i][j-1][l] + ((alpha*dt)/(dr**2))*T[i][j-2][l] + ((-2*alpha*dt)/(dz**2))*T[i][j][l-1] + ((alpha*dt)/(dz**2))*T[i][j][l-2]
            elif j == 0:    #z end and r start
                T[i+1][j][l] = (1 - (alpha*dt)/(r[j]*dr) + (alpha*dt)/(dr**2) + (alpha*dt)/(dz**2))*T[i][j][l] + ((alpha*dt)/(r[j]*dr) - (2*alpha*dt)/(dr**2))*T[i][j+1][l] + ((alpha*dt)/(dr**2))*T[i][j+2][l] + ((-2*alpha*dt)/(dz**2))*T[i][j][l-1] + ((alpha*dt)/(dz**2))*T[i][j][l-2]
            else:           #z end and r middle
                T[i + 1][j][l] = (1 - (2*alpha*dt)/(dr**2) + (alpha*dt)/(dz**2))*T[i][j][l] + ((alpha*dt)/(dr**2) - (alpha*dt)/(2*r[j]*dr))*T[i][j-1][l] + ((alpha*dt)/(dr**2) + (alpha*dt)/(2*r[j]*dr))*T[i][j+1][l] + ((-2*alpha*dt)/(dz**2))*T[i][j][l-1] + ((alpha*dt)/(dz**2))*T[i][j][l-2]
        elif l == 0:
            if j == nr-1:   #z start and r end
                T[i + 1][j][l] = (2 * h_m * dt * Tc) / (rho * cp * dz) + (1 + (k*dt)/(rho*cp*(dr)**2) - (k*dt)/(rho*cp*(dz)**2) - (2*h_m*dt)/(rho*cp*dz))*T[i][j][l] + ((-k*dt)/(rho*cp*(dr)**2))*T[i][j-1][l] + ((k*dt)/(rho*cp*(dz)**2))*T[i][j][l+1]
            elif j == 0:    #z start and r start
                T[i + 1][j][l] = (2*h_m*dt*Tc)/(rho*cp*dz) + (1 - (k*dt)/(rho*cp*(dr)**2) - (k*dt)/(rho*cp*(dz)**2) - (2*h_m*dt)/(rho*cp*dz))*T[i][j][l] + ((k*dt)/(rho*cp*(dr)**2))*T[i][j+1][l] + ((k*dt)/(rho*cp*(dz)**2))*T[i][j][l+1]
            else:           #z start and r middle
                T[i + 1][j][l] = (2 * h_m * dt * Tc) / (rho * cp * dz) + (1 - (k*dt)/(rho*cp*(dz)**2) - (2*h_m*dt)/(rho*cp*dz))*T[i][j][l] + ((k*dt)/(rho*cp*(dr)**2))*T[i][j+1][l] + ((-k*dt)/(rho*cp*(dr)**2))*T[i][j-1][l] + ((k*dt)/(rho*cp*(dz)**2))*T[i][j][l+1]
        else:
            if j == nr-1:   #z middle and r end
                T[i+1][j][l] = (1+(alpha*dt)/(r[j]*dr) + (alpha*dt)/(dr**2) - (2*alpha*dt)/(dz**2))*T[i][j][l] + ((-2*alpha*dt)/(dr**2) - (alpha*dt)/(r[j]*dr))*T[i][j-1][l] + ((alpha*dt)/(dr**2))*T[i][j-2][l] + ((alpha*dt)/(dz**2))*T[i][j][l-1] + ((alpha*dt)/(dz**2))*T[i][j][l+1]

            elif j == 0:    #z middle and r start
                T[i+1][j][l] = (1- (alpha*dt)/(r[j]*dr) + (alpha*dt)/(dr**2)  - (2*alpha*dt)/(dz**2))*T[i][j][l] + ((alpha*dt)/(r[j]*dr) - (2*alpha*dt)/(dr**2))*T[i][j+1][l] + ((alpha*dt)/(dr**2))*T[i][j+2][l] + ((alpha*dt)/(dz**2))*T[i][j][l-1] + ((alpha*dt)/(dz**2))*T[i][j][l+1]
            else:           #z middle and r middle
                T[i + 1][j][l] = (1 - (2*alpha*dt)/(dr**2) - (2*alpha*dt)/(dz**2))*T[i][j][l] + ((alpha*dt)/(dr**2) - (alpha*dt)/(2*r[j]*dr))*T[i][j-1][l] + ((alpha*dt)/(dr**2) + (alpha*dt)/(2*r[j]*dr))*T[i][j+1][l] + ((alpha*dt)/(dz**2))*T[i][j][l-1] + ((alpha*dt)/(dz**2))*T[i][j][l+1]
print(T)
