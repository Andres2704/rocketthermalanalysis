import numpy as np
import matplotlib.pyplot as plt
from constants import *

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

