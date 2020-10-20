import numpy as np
from materials import *
import matplotlib.pyplot as plt
import math
class motor_case():
    def __init__(self, material_case, material_liner,  insulator_thk, case_thk, ri, r_steps, hm, Tc, Ta,
                 burn_time, t_steps, type):
        self.rho_case, self.k_case, self.cp_case = case_selector(material_case)
        self.rho_insulator, self.k_insulator, self.cp_insulator = insulator_selector(material_liner)
        self.insulator_thk = insulator_thk
        self.case_thk = case_thk
        self.r_1 = ri
        self.sectionsr = r_steps
        self.h_m = hm
        self.Tc = Tc
        self.Ta = Ta
        self.t = burn_time
        self.nt = t_steps
        self.alpha_case = float(self.k_case / (self.cp_case * self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))

        if type == 1:
            self.T, self.Values = self.run_analysis_case()
        elif type == 2:
            self.T, self.Values = self.run_analysis_explicit()

    def run_analysis_case(self):
        # Radial coordinate
        self.r_2 = self.r_1 + self.insulator_thk  # Radial position of insulation/casing interface [m]
        self.r_3 = self.r_2 + self.case_thk  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr  # r variation [m]
        self.nr = self.sectionsr + 1  # Nº of radial points

        # Time
        self.dt = self.t / self.nt  # Time variation [s]

        # we have to set the initial and boundary conditions
        self.r = [self.r_1 + i * self.dr for i in range(self.nr)] # Defining the r vector

        A = self.create_Amatrix()

        # Creating the T matrix
        T_implicit = np.zeros([self.nt, self.nr])

        # Setting initial values for the T matrix
        for i in range(self.nr):
            T_implicit[0][i] = self.Ta

        # Getting the inverse matrix of A
        A_inverse = np.linalg.inv(A)

        # Calculating the implicit T matrix
        for i in range(self.nt - 1):
            T_implicit[i][0] = T_implicit[i][0] + (self.h_m*2*self.dt*(self.Tc-T_implicit[i][0]))/(self.rho_insulator*self.cp_insulator*self.dr)
            T_implicit[i + 1] = (np.matmul(A_inverse, T_implicit[i]))

        values = T_implicit[self.nt-1, ::]
        return T_implicit, values

    def run_analysis_explicit(self):
        # Radial coordinate
        self.r_2 = self.r_1 + self.insulator_thk  # Radial position of insulation/casing interface [m]
        self.r_3 = self.r_2 + self.case_thk  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr  # r variation [m]
        self.nr = self.sectionsr + 1  # Nº of radial points

        # Time
        self.dt = self.dr ** 2 / (2 * self.alpha_case)  # Time variation [s]
        self.nt_ex = self.nr
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



        T = np.zeros([self.nt_ex, self.nr])
        for i in range(self.nr):
            T[0][i] = self.Ta

        # Calculating the T-matrix
        for j in range(self.nt_ex - 1):
            for i in range(self.nr):
                if self.r[i] < self.r_2:  # EPDM
                    alpha = self.alpha_insulator
                    rho = self.rho_insulator
                    cp = self.cp_insulator
                    k = self.k_insulator
                if self.r[i] >= self.r_2:  # Corrections are needed here
                    alpha = self.alpha_case
                    rho = self.rho_case
                    cp = self.cp_case
                    k = self.k_case
                if i == self.nr - 1:
                    T[j + 1][i] = T[j][i] + alpha * self.dt * (
                                ((1 / (self.r[i] * self.dr)) * (T[j][i] - T[j][i - 1])) + (
                                    (T[j][i] + T[j][i - 2] - 2 * T[j][i - 1]) /
                                    (self.dr ** 2)))
                elif i == 0:
                    T[j + 1][i] = T[j][i] + (2 / (rho * cp * self.dr)) * (
                                self.h_m * self.dr * (self.Tc - T[j][0]) + k * (T[j][1] - T[j][0]))
                else:
                    T[j + 1][i] = T[j][i] + alpha * self.dt * (
                                ((1 / (2 * self.r[i] * self.dr)) * (T[j][i + 1] - T[j][i - 1])) + (
                                    (T[j][i + 1] + T[j][i - 1] - 2 * T[j][i]) /
                                    (self.dr ** 2)))
        values = T[self.nt_ex - 1, ::]
        return T, values

    def create_Amatrix(self):
        # Creating the A matrix
        A = np.zeros([self.nr, self.nr])

        # Building the A matrix
        A[0][0] = 1 + (2 * self.dt * self.k_insulator) / (self.rho_insulator * self.cp_insulator * self.dr ** 2)
        A[0][1] = -self.k_insulator * 2 * self.dt / (self.dr * self.dr * self.rho_insulator * self.cp_insulator)
        A[self.nr - 1][self.nr - 1] = 1 - self.alpha_case * self.dt / (self.r[self.nr - 1] * self.dr) - self.alpha_case * self.dt / (self.dr ** 2)
        A[self.nr - 1][self.nr - 2] = self.alpha_case * self.dt / (self.r[self.nr - 1] * self.dr) + 2 * self.alpha_case * self.dt / (self.dr ** 2)
        A[self.nr - 1][self.nr - 3] = -self.alpha_case * self.dt / (self.dr ** 2)

        for i in range(1, self.nr - 1):
            if self.r[i] >= self.r_2:
                alpha = self.alpha_case
            else:
                alpha = self.alpha_insulator
            A[i][i - 1] = -(alpha * self.dt) / (self.dr ** 2) + alpha * self.dt / (2 * self.dr * self.r[i])
            A[i][i] = 1 + 2 * alpha * self.dt / (self.dr ** 2)
            A[i][i + 1] = -(alpha * self.dt) / self.dr ** 2 - alpha * self.dt / (2 * self.dr * self.r[i])
        return A

class bulkhead(object):
    def __init__(self, material_case, material_liner, insulator_thk, bulkhead_thk, radius_inner, radius_outer, burn_time, sectionsr, hm, Ta, Tc ):
        # Explicit atributes for bulkhead
        self.rho_case, self.k_case, self.cp_case = case_selector(material_case) # Case properties
        self.rho_insulator, self.k_insulator, self.cp_insulator = insulator_selector(material_liner) # Insulator properties
        self.h_m = hm
        self.insulator_thk = insulator_thk
        self.z = bulkhead_thk
        self.radius_inner = radius_inner
        self.radius_outer = radius_outer
        self.t = burn_time
        self.sectionsr = sectionsr
        self.Ta = Ta
        self.Tc = Tc
        self.alpha_case = float(self.k_case/(self.cp_case*self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))

        self.T_bulkhead = self.run_analysis_bulkhead()

    def run_analysis_bulkhead(self):

        nr = self.sectionsr + 1
        dr = (self.radius_outer) / self.sectionsr  # r variation [m]
        r_1 = dr  # Radial position of insulation beginning [m]
        r_2 = r_1 + self.radius_inner  # Radial position of insulation/casing interface [m]
        r_3 = self.radius_outer  # Radial position of casing end [m]

        z_1 = 0
        z_2 = self.insulator_thk
        z_3 = z_2 + self.z
        dz = dr / 1.4  # z variation [m]
        nz = math.ceil(self.z / dz)  # Nº of z points

        dt = ((dr ** 2) * (dz ** 2)) / (1.99 * self.alpha_case * (dr ** 2 + dz ** 2))
        nt = math.ceil(self.t / dt)

        d_insulator = self.alpha_insulator * (dt / (dr ** 2))  # Insulation Fourier number
        d_case = self.alpha_case * (dt / (dr ** 2))  # Casing Fourier number

        # we have to set the initial/boundary contitions
        r = [r_1 + i * dr for i in range(nr)]
        z = [z_1 + i * dr for i in range(nz)]

        T = np.zeros([nt, nr, nz])
        for i in range(nr):
            for j in range(nz):
                T[0][i][j] = self.Ta

        # Calculating the T-matrix
        for i in range(nt - 1):
            for l in range(nz):
                for j in range(nr):
                    if z[l] <= z_2:  # Set insulation material
                        alpha = self.alpha_insulator
                        rho = self.rho_insulator
                        cp = self.cp_insulator
                        k = self.k_insulator
                    if z[l] > z_2:  # Set casing material
                        alpha = self.alpha_case
                        rho = self.rho_case
                        cp = self.cp_case
                        k = self.k_case
                    if l == nz - 1:
                        if j == nr - 1:  # z end and r end
                            T[i + 1][j][l] = (1 + (alpha * dt) / (r[j] * dr) + (alpha * dt) / (dr ** 2) + (
                                        alpha * dt) / (dz ** 2)) * T[i][j][l] + (
                                                         (-2 * alpha * dt) / (dr ** 2) - (alpha * dt) / (r[j] * dr)) * \
                                             T[i][j - 1][l] + ((alpha * dt) / (dr ** 2)) * T[i][j - 2][l] + (
                                                         (-2 * alpha * dt) / (dz ** 2)) * T[i][j][l - 1] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l - 2]
                        elif j == 0:  # z end and r start
                            T[i + 1][j][l] = (1 - (alpha * dt) / (r[j] * dr) + (alpha * dt) / (dr ** 2) + (
                                        alpha * dt) / (dz ** 2)) * T[i][j][l] + (
                                                         (alpha * dt) / (r[j] * dr) - (2 * alpha * dt) / (dr ** 2)) * \
                                             T[i][j + 1][l] + ((alpha * dt) / (dr ** 2)) * T[i][j + 2][l] + (
                                                         (-2 * alpha * dt) / (dz ** 2)) * T[i][j][l - 1] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l - 2]
                        else:  # z end and r middle
                            T[i + 1][j][l] = (1 - (2 * alpha * dt) / (dr ** 2) + (alpha * dt) / (dz ** 2)) * T[i][j][
                                l] + ((alpha * dt) / (dr ** 2) - (alpha * dt) / (2 * r[j] * dr)) * T[i][j - 1][l] + (
                                                         (alpha * dt) / (dr ** 2) + (alpha * dt) / (2 * r[j] * dr)) * \
                                             T[i][j + 1][l] + ((-2 * alpha * dt) / (dz ** 2)) * T[i][j][l - 1] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l - 2]
                    elif l == 0:
                        if j == nr - 1:  # z start and r end
                            T[i + 1][j][l] = (2 * self.h_m * dt * Tc) / (rho * cp * dz) + (
                                        1 + (k * dt) / (rho * cp * (dr) ** 2) - (k * dt) / (rho * cp * (dz) ** 2) - (
                                            2 * self.h_m * dt) / (rho * cp * dz)) * T[i][j][l] + (
                                                         (-k * dt) / (rho * cp * (dr) ** 2)) * T[i][j - 1][l] + (
                                                         (k * dt) / (rho * cp * (dz) ** 2)) * T[i][j][l + 1]
                        elif j == 0:  # z start and r start
                            T[i + 1][j][l] = (2 * self.h_m * dt * Tc) / (rho * cp * dz) + (
                                        1 - (k * dt) / (rho * cp * (dr) ** 2) - (k * dt) / (rho * cp * (dz) ** 2) - (
                                            2 * self.h_m * dt) / (rho * cp * dz)) * T[i][j][l] + (
                                                         (k * dt) / (rho * cp * (dr) ** 2)) * T[i][j + 1][l] + (
                                                         (k * dt) / (rho * cp * (dz) ** 2)) * T[i][j][l + 1]
                        else:  # z start and r middle
                            T[i + 1][j][l] = (2 * self.h_m * dt * Tc) / (rho * cp * dz) + (
                                        1 - (k * dt) / (rho * cp * (dz) ** 2) - (2 * self.h_m * dt) / (rho * cp * dz)) * \
                                             T[i][j][l] + ((k * dt) / (rho * cp * (dr) ** 2)) * T[i][j + 1][l] + (
                                                         (-k * dt) / (rho * cp * (dr) ** 2)) * T[i][j - 1][l] + (
                                                         (k * dt) / (rho * cp * (dz) ** 2)) * T[i][j][l + 1]
                    else:
                        if j == nr - 1:  # z middle and r end
                            T[i + 1][j][l] = (1 + (alpha * dt) / (r[j] * dr) + (alpha * dt) / (dr ** 2) - (
                                        2 * alpha * dt) / (dz ** 2)) * T[i][j][l] + (
                                                         (-2 * alpha * dt) / (dr ** 2) - (alpha * dt) / (r[j] * dr)) * \
                                             T[i][j - 1][l] + ((alpha * dt) / (dr ** 2)) * T[i][j - 2][l] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l - 1] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l + 1]

                        elif j == 0:  # z middle and r start
                            T[i + 1][j][l] = (1 - (alpha * dt) / (r[j] * dr) + (alpha * dt) / (dr ** 2) - (
                                        2 * alpha * dt) / (dz ** 2)) * T[i][j][l] + (
                                                         (alpha * dt) / (r[j] * dr) - (2 * alpha * dt) / (dr ** 2)) * \
                                             T[i][j + 1][l] + ((alpha * dt) / (dr ** 2)) * T[i][j + 2][l] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l - 1] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l + 1]
                        else:  # z middle and r middle
                            T[i + 1][j][l] = (1 - (2 * alpha * dt) / (dr ** 2) - (2 * alpha * dt) / (dz ** 2)) * \
                                             T[i][j][l] + ((alpha * dt) / (dr ** 2) - (alpha * dt) / (2 * r[j] * dr)) * \
                                             T[i][j - 1][l] + (
                                                         (alpha * dt) / (dr ** 2) + (alpha * dt) / (2 * r[j] * dr)) * \
                                             T[i][j + 1][l] + ((alpha * dt) / (dz ** 2)) * T[i][j][l - 1] + (
                                                         (alpha * dt) / (dz ** 2)) * T[i][j][l + 1]

        return T

if __name__ == "__main__":
    # CASING ANALYSIS
    material_case = 'Aluminum 6061 - T6'  # Case material ['aluminium6061t6']
    material_liner = 'EPDM'     # Insulator Material ['epdm']
    ri = 0.11                   # Initial cylinder radius
    insulatior_thk = 0.003      # Insulator thickness
    case_thk = 0.00555          # Case thickness
    r_steps = 100               # Radial steps
    hm = 1295                   # Average convection coefficient
    Tc = 1600                   # Combustion chamber temperature
    Ta = 297                    # Ambient temperature
    burn_time = 7               # Burn time
    t_steps = 500               # time steps

    """
    type = 2
    motor1 = motor_case(material_case, material_liner, insulatior_thk, case_thk, ri, r_steps, hm, Tc, Ta, burn_time, t_steps, type)

    print(motor1.T)
    
    # Building graphs
    r_a = np.array(motor.r)
    locals = [r_a >= motor.r_2]
    r_a = r_a[tuple(locals)]

    fig, ax = plt.subplots(figsize=(10,25))
    ax.plot(r_a, motor.Values[motor.r.index(r_a[0])::], 'r-+', linewidth=2, label='Implicit Method')
    ax.plot(r_a, motor1.Values[motor.r.index(r_a[0])::], 'b', linewidth=2, label='Explicit Method')
    ax.set_title('Temperature distribuition along the case - Insulation tichkness: ' + str(insulatior_thk*1000) + ' mm')
    ax.set_xlabel('Radial position [m]')
    ax.set_ylabel('Temperature [K]')
    ax.legend()
    plt.show()
    """
    
    # BULKHEAD ANALYSIS
    radius_inner = 0.064
    radius_outer = 0.07
    bulkhead_thk = 0.007

    bulk = bulkhead(material_case, material_liner, insulatior_thk, bulkhead_thk, radius_inner, radius_outer, burn_time, r_steps, hm, Ta, Tc)

    print(bulk.T_bulkhead)

