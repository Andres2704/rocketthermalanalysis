import numpy as np
np.set_printoptions(precision=2)
# Liner properties
k_insulator = 0.2                                                      # Insulation thermal conductivity [W/m-K]
rho_insulator = 860                                                    # Insulation density [kg/m3]
cp_insulator = 2000                                                    # Insulation specific heat capacity [J/kg-K]
alpha_insulator = float(k_insulator / (cp_insulator * rho_insulator))  # Insulation thermal diffusivity [m^2/s]
# ------------------------------------------------
# Case properties
k_case = 167                                                           # Casing thermal conductivity [W/m-K]
rho_case = 2700                                                        # Casing density [kg/m3]
cp_case = 896                                                          # Casing specific heat capacity [J/kg-K]
alpha_case = float(k_case / (cp_case * rho_case))                      # Casing thermal diffusivity [m^2/s]

hm = 1295


# ------------------------------------------------
# Time Domain
t_b = 5.5                                                              # Burn time
Nt = 500                                                               # Number of time steps
dt = t_b/Nt                                                            # Time steps
# ------------------------------------------------
# Radial domain
liner_thk = 0.003                                                      # Insulation thickness [m]
case_thk = 0.00538                                                     # Case thickness |VERIFY THIS VALUEEEE| [m]
r1 = 0.11                                                              # Radial position of insulation beginning [m]
r2 = r1 + liner_thk                                                    # Radial position of insulation/casing interface [m]
r3 = r2 + case_thk                                                     # Radial position of casing end [m]
Nr = 50                                                                # Number of radial steps
dr = (r3-r1)/Nr                                                        # Radial steps
# ------------------------------------------------
# Fourier parameters
d_insulator = alpha_insulator * (dt / (dr*dr))                       # Insulation Fourier number
d_case =  alpha_case * (dt / (dr*dr))                                 # Case Fourier number
co_case = alpha_case*(dt/dr)
co_liner = alpha_insulator*(dt/dr)

# ------------------------------------------------
# Adding some constants an defining the r vector due the discretization and the T-matrix
Tc = 1600                                                              # Temperature at the combustion chamber [K]
Ta = 297                                                               # Average ambient temperature at Tatúi-Sp.
T = np.zeros([Nr+1, Nt+1])

r = []
t = []

# ------------------------------------------------
# Setting the initial conditions
for n in range(Nr+1):
    T[n, 0] = Ta                                                      # The initial temperature of the liner at t=0
    r.append(r1 + n*dr)


for j in range(Nt+1):
    t.append(j*dt)

# Setting the boundary conditions
for i in range(Nr):
    for j in range(Nt):
        if i==Nr:
            T[i, j+1] = T[i,j] + alpha_case*dt*((1/(r[i]*dr)) + (T[i,j]-T[i-1,j]) + ((T[i,j] + T[i-2, j] - 2*T[i-1, j])/dr*dr))
        if i==0:
            T[i, j+1] = T[i, j] + (2 / (rho_insulator * cp_insulator * dr)) * (hm * dt * (Tc - T[0, j]) + k_insulator * (T[1, j] - T[0, j]))

print(T.transpose())
# ------------------------------------------------
# Defining the left(MMl) and right(MMr) matrices of Crank-Nicholson method
# The MMr are tridiagonal matrices, so for that we will define the three diagonal first
aal = [-d_insulator*(r[r_i]/(r[r_i]-co_liner)) for r_i in range(0,Nr-2) if r[r_i]<=r2]                              # Below the main diagonal MMl
aal = aal + [-d_case*(r[r_i]/(r[r_i]-co_case)) for r_i in range(0,Nr-2) if r[r_i]>r2]                               # Below the main diagonal MMl
bbl = [(2+2*d_insulator + co_liner/r[r_i])*(r[r_i]/(r[r_i]-co_liner)) for r_i in range(0,Nr-1) if r[r_i]<=r2]       # Main diagonal in MMl
bbl = bbl + [(2+2*d_case + co_case/r[r_i])*(r[r_i]/(r[r_i]-co_case)) for r_i in range(0, Nr-1) if r[r_i]>r2]        # Main diagonal in MMl
ccl = [-d_insulator*(r[r_i]/(r[r_i]-co_liner)) for r_i in range(0,Nr-2) if r[r_i]<=r2]                              # Above the main matrix MMl
ccl = ccl + [-d_case*(r[r_i]/(r[r_i]-co_case)) for r_i in range(0,Nr-2) if r[r_i]>r2]                               # Above the main matrix MMl
MMl = np.diagflat(bbl,0) + np.diagflat(aal,-1) + np.diagflat(ccl,1)

aar = [d_insulator*(r[r_i]/(r[r_i]-co_liner)) for r_i in range(0,Nr-2) if r[r_i]<=r2]                               # Below the main diagonal MMr
aar = aar + [d_case*(r[r_i]/(r[r_i]-co_case)) for r_i in range(0,Nr-2) if r[r_i]>r2]                                # Below the main diagonal MMr
bbr = [(2-2*d_insulator - co_liner/r[r_i])*(r[r_i]/(r[r_i]-co_liner)) for r_i in range(0,Nr-1) if r[r_i]<=r2]       # Main diagonal in MMr
bbr = bbr + [(2-2*d_case - co_case/r[r_i])*(r[r_i]/(r[r_i]-co_case)) for r_i in range(0, Nr-1) if r[r_i]>r2]        # Main diagonal in MMr
ccr = [d_insulator*(r[r_i]/(r[r_i]-co_liner)) for r_i in range(0,Nr-2) if r[r_i]<=r2]                               # Above the main matrix MMr
ccr = ccr + [d_case*(r[r_i]/(r[r_i]-co_case)) for r_i in range(0,Nr-2) if r[r_i]>r2]                                # Above the main matrix MMr
MMr = np.diagflat(bbr,0) + np.diagflat(aar,-1) + np.diagflat(ccr,1)

#Applying the Crank-Nicholson method
for j in range(Nt):
    A = np.linalg.inv(MMl)
    B = np.matmul(A, MMr)
    T[1:Nr, j+1] = np.matmul(B, T[1:Nr, j])


"""
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 10))
plt.plot(r, T[:,0], label='T=0')
plt.plot(r, T[:,round(Nt/100)], label='T=Nt/100')
plt.plot(r, T[:,round(Nt/10)], label='T=Nt/10')
plt.plot(r, T[:,round(Nt)], label='T=Nt')
plt.show()


import numpy as np
cond =1/2
L = 1
T = 1

Nt = 2500
Dt = T/Nt
Nr = 50
Dr = L/Nr
b = cond*Dt/(Dr*Dr)
r = []

u = np.zeros([Nr+1, Nt+1])
for n in range(0,Nr+1):
    r.append((n*Dr))
    u[n, 0] = np.sin(np.pi*r[n])


t = []
for j in range(0, Nt+1):
    u[0, j] = 0
    u[Nr, j] = 0
    t.append(j*Dt)



aal = [-b for i in range(0, Nr-2)]         # Below the main diagonal MMl
bbl = [2+2*b for i in range(0,Nr-1)]      # Main diagonal in MMl
ccl = [-b for i in range(0, Nr-2)]       # Above the main matrix MMl
MMl = np.diagflat(bbl,0) + np.diagflat(aal,-1) + np.diagflat(ccl,1)

aar = [b for i in range(0, Nr-2)]         # Below the main diagonal MMl
bbr = [2+2*(-b) for i in range(0,Nr-1)]      # Main diagonal in MMl
ccr = [b for i in range(0, Nr-2)]       # Above the main matrix MMl
MMr = np.diagflat(bbr,0) + np.diagflat(aar,-1) + np.diagflat(ccr,1)



for j in range(0, Nt):
    A = np.linalg.inv(MMl)
    B = np.matmul(A, MMr)
    u[1:Nr, j + 1] = np.matmul(B, u[1:Nr, j])


import matplotlib.pyplot as plt

plt.plot(r, u[:,0],  'b', label='t=0')
plt.plot(r, u[:,round(Nt/100)],  'g--', label='t=Nt/100')
plt.plot(r, u[:,round(Nt/10)], 'r', label='t=Nt/10')
plt.plot(r, u[:,Nt], 'yo', label='t=Nt')
plt.legend()
plt.show()
"""
"""
########################### - Convergence analysis
import numpy as np
import matplotlib.pyplot as plt
from constants import *

# External for with the porpouse to analyse the convergence of the methos
steps = [i for i in range(10, 2000, 10)]
conv_i, conv = [], []

for typ in range(2):
    for step_i in steps:
      if typ == 0:
          nt = step_i
          dt = t / nt
      if typ == 1:
          sectionsr = step_i  # Nº of sections between r_1 and r_3
          dr = (r_3 - r_1) / sectionsr  # r variation [m]
          nr = sectionsr + 1  # Nº of radial points

      Fo_insulator = alpha_insulator * (dt / (dr ** 2))                                           # Insulation Fourier number
      Fo_case = alpha_case * (dt / (dr ** 2))                                                     # Casing Fourier number

      Bi_insulator = h_m * (dr) / k_insulator                                                     # Insulation Biot number
      Bi_case = h_m * (dr) / k_case                                                               # Casing Biot number

      # we have to set the initial/boundary contitions
      r = [r_1 + i * dr for i in range(nr)]
      z = [i * dz for i in range(nz)]

      #Creating the A matrix
      A = np.zeros([nr, nr])

      #Building the A matrix
      A[0][0] = 1 + (2*dt*k_insulator)/(rho_insulator*cp_insulator*dr**2)
      A[0][1] = -k_insulator*2*dt / (dr*dr*rho_insulator*cp_insulator)
      A[nr-1][nr-1] = 1 - alpha_case * dt / (r[nr-1] * dr) - alpha_case * dt / (dr ** 2)
      A[nr-1][nr-2] = alpha_case * dt / (r[nr-1] * dr) + 2 * alpha_case * dt / (dr ** 2)
      A[nr-1][nr-3] = -alpha_case * dt / (dr ** 2)

      for i in range(1, nr-1):
        if r[i] >= r_2: alpha = alpha_case
        else: alpha = alpha_insulator
        A[i][i - 1] = -(alpha * dt) / (dr ** 2) + alpha * dt / (2 * dr * r[i])
        A[i][i] = 1 + 2 * alpha * dt / (dr ** 2)
        A[i][i + 1] = -(alpha * dt) / dr ** 2 - alpha * dt / (2 * dr * r[i])

      #Creating the T matrix
      T = np.zeros([nt, nr])


      #Setting initial values for the T matrix
      for i in range(nr):
          T[0][i] = Ta

      #Getting the inverse matrix of A
      A_inverse = np.linalg.inv(A)


      # Calculating the T matrix
      for i in range(nt-1):
          T[i][0] = T[i][0] + 2*dt*h_m*(Tc-T[i][0])/(rho_insulator*cp_insulator*dr)
          T[i+1] = (np.matmul(A_inverse, T[i]))

      conv_i.append(T[nt-1,nr-1])

    conv.append(conv_i)
    conv_i = []

fig, ax = plt.subplots(1,2, figsize=(12,10))
ax[0].plot(steps, conv[0])
ax[0].set_xlabel('Time steps')
ax[0].set_ylabel('Final temperature at the case [K]')
ax[1].plot(steps, conv[1])
ax[1].set_xlabel('Radial steps')
ax[1].set_ylabel('Final temperature at the case [K]')
fig2, ax1 = plt.subplots(figsize=(7,15))
ax1.plot(r[9:], T[99,9:nr], label='Burn time: 7s')
ax1.set_title('Temperature distribution at motor case along the radial direction')
ax1.legend()
plt.show()

#########################
Temperature graphs
r_a = np.array(r)
locals = [r_a >= r_2]
r_a = r_a[tuple(locals)]

fig2, ax1 = plt.subplots(figsize=(7,15))
ax1.plot(r_a, T[nt-1,r.index(r_a[0]):nr], 'r', label='Burn time: '+str(t)+'s')
ax1.set_title('Temperature distribution at motor case along the radial direction')
ax1.set_xlabel('Radial position of motor case [m]')
ax1.set_ylabel('Temperature [K]')
ax1.legend()
plt.show()

print(T) 



###########################

 def run_analysis_explicit(self):
        # Radial coordinate
        self.r_1 = 0.11  # Radial position of insulation beginning [m]
        self.r_2 = self.r_1 + self.insulator_thk  # Radial position of insulation/casing interface [m]
        self.r_3 = self.r_2 + self.case_thk  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr  # r variation [m]
        self.nr = self.sectionsr + 1  # Nº of radial points

        # Time
        self.dt = self.t / self.nt  # Time variation [s]

        # we have to set the initial and boundary conditions
        self.r = [self.r_1 + i * self.dr for i in range(self.nr)]  # Defining the r vector

        d_insulator = self.alpha_insulator * (self.dt / (self.dr ** 2))  # Insulation Fourier number
        d_case = self.alpha_case * (self.dt / (self.dr ** 2))  # Casing Fourier number


        if d_insulator > 0.5 or d_case > 0.5:
            print('Stability criteria was not achieved, the following values must be less than 0.5')
            print('Fourier number case: ', d_case)
            print('Fourier number insulator: ', d_insulator)
            print('Try to decrease the radial steps')
            return 0

        # Height coordinate
        z = 1  # Height of combustion chamber [m]
        sectionsz = 1  # Nº of sections beewtween 0 and z
        self.dz = z / sectionsz  # z variation [m]
        nz = sectionsz + 2  # Nº of z points

        T = np.zeros([self.nt, self.nr, nz])
        for i in range(self.nr):
            for j in range(nz):
                T[0][i][j] = Ta

        # Calculating the T-matrix
        for j in range(self.nt - 1):
            for i in range(self.nr):
                for l in range(nz):
                    if self.r[i] <= self.r_2:  # Set insulation material
                        alpha = self.alpha_insulator
                        rho = self.rho_insulator
                        cp = self.cp_insulator
                        k = self.k_insulator
                    if self.r[i] > self.r_2:  # Set casing material
                        alpha = self.alpha_case
                        rho = self.rho_case
                        cp = self.cp_case
                        k = self.k_case
                    if i == self.nr - 1:
                        if l == nz - 1:  # r end and z end
                            T[j + 1][i][l] = T[j][i][l] + alpha * self.dt * (((T[j][i][l] - 2 * T[j][i][l - 1] + T[j][i][l - 2]) / (self.dz ** 2)) +
                                            ((T[j][i][l] - T[j][i - 1][l]) / (self.r[i] * self.dr))+((T[j][i][l] + T[j][i - 2][l] - 2 * T[j][i - 1][l]) / (self.dr ** 2)))
                        elif l == 0:  # r end and z start
                            T[j + 1][i][l] = T[j][i][l] + alpha * self.dt * (((T[j][i][l] - 2 * T[j][i][l + 1] + T[j][i][l + 2]) / (self.dz ** 2)) + (
                                            (T[j][i][l] - T[j][i - 1][l]) / (self.r[i] * self.dr)) + ((T[j][i][l] + T[j][i - 2][l] - 2 * T[j][i - 1][l]) / (self.dr ** 2)))
                        else:  # r end and z middle
                            T[j + 1][i][l] = T[j][i][l] + alpha * self.dt * (((T[j][i][l + 1] - 2 * T[j][i][l] + T[j][i][l - 1]) / (self.dz ** 2)) + (
                                            (T[j][i][l] - T[j][i - 1][l]) / (self.r[i] * self.dr)) + ((T[j][i][l] + T[j][i - 2][l] - 2 * T[j][i - 1][l]) / (self.dr ** 2)))
                    elif i == 0:
                        if l == nz - 1:  # r start and z end
                            T[j + 1][i][l] = T[j][i][l] + ((2 * self.dt) / (rho * cp * self.dr * self.dz)) * (self.h_m * self.dz * (Tc - T[j][i][l]) + k * self.dz * (
                                            T[j][i + 1][l] - T[j][i][l]) / self.dr + k * self.dr * (T[j][i][l] - T[j][i][l - 1]) / (2 * self.dz))
                        elif l == 0:  # r start and z start
                            T[j + 1][i][l] = T[j][i][l] + ((2 * self.dt) / (rho * cp * self.dr * self.dz)) * (self.h_m * self.dz * (Tc - T[j][i][l]) + k * self.dz * (
                                            T[j][i + 1][l] - T[j][i][l]) / self.dr + k * self.dr * (T[j][i][l + 1] - T[j][i][l]) / (2 * self.dz))
                        else:  # r start and z middle
                            T[j + 1][i][l] = T[j][i][l] + ((2 * self.dt) / (rho * cp * self.dr * self.dz)) * (
                                        self.h_m * self.dz * (Tc - T[j][i][l]) + k * self.dz * (T[j][i + 1][l] - T[j][i][l]) / self.dr + k * self.dr *
                                        (T[j][i][l + 1] - T[j][i][l - 1]) / (2 * self.dz))
                    else:
                        if l == nz - 1:  # r middle and z end
                            T[j + 1][i][l] = T[j][i][l] + alpha * self.dt * (((T[j][i + 1][l] - T[j][i - 1][l]) / (2 * self.r[i] * self.dr)) + ((T[j][i + 1][l] +
                                                                               T[j][i - 1][l - 2] - 2 * T[j][i][l]) / (self.dr ** 2)) +
                                                                             ((T[j][i][l] + T[j][i][l - 2] - 2 * T[j][i][l - 1]) / (self.dz ** 2)))
                        elif l == 0:  # r middle and z start
                            T[j + 1][i][l] = T[j][i][l] + alpha * self.dt * (((T[j][i + 1][l] - T[j][i - 1][l]) / (2 * self.r[i] * self.dr)) + ((T[j][i + 1][l] +
                                                                              T[j][i - 1][l - 2] - 2 * T[j][i][l]) / (self.dr ** 2)) +
                                                                             ((T[j][i][l] + T[j][i][l + 2] - 2 * T[j][i][l + 1]) / (self.dz ** 2)))
                        else:  # r middle and z middle
                            T[j + 1][i][l] = T[j][i][l] + alpha * self.dt * (((T[j][i + 1][l] - T[j][i - 1][l]) / (2 * self.r[i] * self.dr)) + ((T[j][i + 1][l] +
                                                                            T[j][i - 1][l - 2] - 2 * T[j][i][l]) / (self.dr ** 2)) + (
                                                                            (T[j][i][l + 1] + T[j][i][l - 1] - 2 * T[j][i][l]) / (self.dz ** 2)))
        return T

"""
