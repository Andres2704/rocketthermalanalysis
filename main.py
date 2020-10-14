import numpy as np
from constants import *

# Asking to the user what kind of analisys it will be running
print('1 - Motor case\n2 - Bulkhead\n')
choose = int(input('Choose: '))

if choose==1:
    # Setting the radial and coordinates, also setting the burn time and steps
    # ------------------------------------------------
    # Inputs
    # Radial coordinate
    insulator_thk = 0.003                                                                       # Insulation thickness [m]
    case_thk = 0.00555                                                                          # Case thickness |VERIFY THIS VALUEEEE| [m]
    r_1 = 0.11                                                                                  # Radial position of insulation beginning [m]
    r_2 = r_1 + insulator_thk                                                                   # Radial position of insulation/casing interface [m]
    r_3 = r_2 + case_thk                                                                        # Radial position of casing end [m]
    sectionsr = int(input('Radial and z steps: '))                                              # Nº of sections between r_1 and r_3
    dr = (r_3 - r_1) / sectionsr                                                                # r variation [m]
    nr = sectionsr + 1                                                                          # Nº of radial points

    # Height coordinate
    z = float(input('Motor lenght [m]: '))                                                      # Height of combustion chamber [m]
    dz = dr                                                                                     # z variation [m]
    sectionsz = math.ceil(z / dz)                                                               # Nº of sections beewtween 0 and z
    nz = sectionsz + 1                                                                          # Nº of z points

    # Time
    t = float(input('Burn time [s]: '))                                                         # Burn time [s]
    nt = int(input('Time steps: '))                                                             # Time steps
    dt = t / nt                                                                                 # Time variation [s]

    # we have to set the initial and boundary contitions
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
