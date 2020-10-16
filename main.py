import numpy as np
from materials import *

class motor_case():
    def __init__(self, material_case, material_liner,  insulator_thk, case_thk, r_steps, hm, Tc, Ta, motor_lenght,
                 burn_time, t_steps, type):
        self.rho_case, self.k_case, self.cp_case = case_selector(material_case)
        self.rho_insulator, self.k_insulator, self.cp_insulator = insulator_selector(material_liner)
        self.insulator_thk = insulator_thk
        self.case_thk = case_thk
        self.sectionsr = r_steps
        self.z = motor_lenght
        self.h_m = hm
        self.Tc = Tc
        self.Ta = Ta
        self.t = burn_time
        self.nt = t_steps
        self.alpha_case = float(self.k_case / (self.cp_case * self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))

        if type == 1:
            self.T = self.run_analysis_case()
        elif type == 2:
            self.T = self.run_analysis_explicit()

    def run_analysis_case(self):
        # Radial coordinate
        self.r_1 = 0.11  # Radial position of insulation beginning [m]
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

        return T_implicit

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

        T = np.zeros([self.nt, self.nr])
        for i in range( self.nr):
            T[0][i] =  self.Ta

        # Calculating the T-matrix
        for j in range( self.nt - 1):
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
                    T[j + 1][i] = T[j][i] + alpha * self.dt * (((1 / (self.r[i] * self.dr)) * (T[j][i] - T[j][i - 1])) + ((T[j][i] + T[j][i - 2] - 2 * T[j][i - 1]) /
                                                                (self.dr ** 2)))
                elif i == 0:
                    T[j + 1][i] = T[j][i] + (2 / (rho * cp * self.dr)) * (self.h_m * self.dr * (Tc - T[j][0]) + k * (T[j][1] - T[j][0]))
                else:
                    T[j + 1][i] = T[j][i] + alpha * self.dt * (((1 / (2 * self.r[i] * self.dr)) * (T[j][i + 1] - T[j][i - 1])) + ((T[j][i + 1] + T[j][i - 1] - 2 * T[j][i]) /
                                                                                                                                  (self.dr ** 2)))
        return T

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
    def __init__(self, material_case, material_liner, insulator_thk, case_thk, r_steps, hm, Tc, Ta, burn_time, t_steps,
                 z_steps, z_thk, r_ex, r_in_ex,):
        # Explicit atributes for bulkhead
        self.rho_case, self.k_case, self.cp_case = case_selector(material_case) # Case properties
        self.rho_insulator, self.k_insulator, self.cp_insulator = insulator_selector(material_liner) # Insulator properties
        self.insulator_thk_ex = insulator_thk  # Insulator thickness
        self.t = burn_time
        self.nt_ex = t_steps    # Time steps
        self.nz = z_steps # Z steps
        self.sectionsr_ex = r_steps # Radial steps
        self.z_thk_ex = z_thk # Z tichnes
        self.r_ex_ex = r_ex
        self.r_in_ex = r_in_ex
        self.alpha_case = float(self.k_case / (self.cp_case * self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))
        self.Tc = Tc # Combustion chamber temperature
        self.Ta = Ta # Ambient temperature
        self.h_m = hm # Average Convective heat coefficient
        self.T = self.run_analysis_bulkhead()

    def run_analysis_bulkhead(self):
        # Radial coordinate
        self.r_1 = 0  # Radial start position [m]
        self.r_2 = self.r_in_ex  # Radial position of insulation end [m]
        self.r_3 = self.r_2 + self.r_ex_ex  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr_ex  # r variation [m]
        self.nr_ex = self.sectionsr_ex + 1  # Nº of radial points

        # Height coordinate
        self.z_1 = 0  #Z start position [m]
        self.z_2 = self.insulator_thk_ex  #Bulkhead and insulator interface position [m]
        self.z_3 = self.z_thk_ex + self.z_2  #Bulkhead end [m]
        self.dz = self.z_3/self.nz  #z variation [m]
        self.nz_ex = self.nz + 1  # Nº of z points

        # Time
        self.dt = self.t / self.nt_ex  # Time variation [s]

        # we have to set the initial and boundary contitions
        self.r = [self.r_1 + i * self.dr for i in range(self.nr_ex)]
        self.z = [self.z_1 + i * self.dz for i in range(self.nz_ex)]


        T_explicit = np.zeros([self.nt_ex, self.nr_ex, self.nz_ex])
        for i in range(self.nr_ex):
            for j in range(self.nz_ex):
                T_explicit[0][i][j] = Ta

        # Calculating the explicit T-matrix
        for j in range(self.nt_ex - 1):
            for i in range(self.nr_ex):
                for l in range(self.nz_ex-1):
                    if self.z[i] <= self.z_2:  # Set insulation material
                        alpha = self.alpha_insulator
                        rho = self.rho_insulator
                        cp = self.cp_insulator
                        k = self.k_insulator
                    if self.z[i] > self.z_2:  # Set casing material
                        alpha = self.alpha_case
                        rho = self.rho_case
                        cp = self.cp_case
                        k = self.k_case
                    if i == self.nr_ex - 1:
                        if l == self.nz_ex - 1:  # r end and z end
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + \
                                                      alpha * self.dt * (((T_explicit[j][i - 1][l] - T_explicit[j][i][l]) / (self.r[i] * self.dr)) +
                                                                         ((T_explicit[j][i - 2][l] + T_explicit[j][i][l] - 2 * T_explicit[j][i - 1][l]) / (self.dr ** 2)) +
                                                                         ((T_explicit[j][i][l - 2] + T_explicit[j][i][l] - 2 * T_explicit[j][i][l - 1]) / (self.dz ** 2)))
                        elif l == 0:  # r end and z start
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + ((2 * self.dt) / (rho * cp * self.dr * self.dz)) * \
                                                      ((k * self.dz * (T_explicit[j][i][l] - T_explicit[j][i - 1][l]) / (2 * self.dr)) +
                                                       (k * self.dr * (T_explicit[j][i][l + 1] - T_explicit[j][i][l]) / self.dz) +
                                                       self.h_m * self.dr * (self.Tc - T_explicit[j][i][0]))
                        else:  # r end and z middle
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + \
                                                      alpha * self.dt * (((T_explicit[j][i][l] - T_explicit[j][i-1][l]) / (self.r[i] * self.dr)) +
                                                                         ((T_explicit[j][i - 2][l] + T_explicit[j][i][l] - 2 * T_explicit[j][i - 1][l]) / (self.dr ** 2)) +
                                                                         ((T_explicit[j][i][l + 1] + T_explicit[j][i][l - 1] - 2 * T_explicit[j][i][l]) / (self.dz ** 2)))
                    elif i == 0:
                        if l == self.nz - 1:  # r start and z end
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + \
                                                      alpha * self.dt * (((T_explicit[j][i + 1][l] - T_explicit[j][i][l]) / (self.r[i] * self.dr)) +
                                                                         ((T_explicit[j][i + 2][l] + T_explicit[j][i][l] - 2 * T_explicit[j][i + 1][l]) / (self.dr ** 2)) +
                                                                         ((T_explicit[j][i][l] + T_explicit[j][i][l - 2] - 2 * T_explicit[j][i][l-1]) / (self.dz ** 2)))
                        elif l == 0:  # r start and z start
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + ((2 * self.dt) / (rho * cp * self.dr * self.dz)) * \
                                                      ((k * self.dz * (T_explicit[j][i + 1][l] - T_explicit[j][i][l]) / (2 * self.dr)) +
                                                       (k * self.dr * (T_explicit[j][i][l + 1] - T_explicit[j][i][l]) / self.dz) +
                                                       self.h_m * self.dr * (self.Tc - T_explicit[j][i][0]))
                        else:  # r start and z middle
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + \
                                                      alpha * self.dt * (((T_explicit[j][i + 1][l] - T_explicit[j][i][l]) / (self.r[i] * self.dr)) +
                                                                         ((T_explicit[j][i + 2][l] + T_explicit[j][i][l] - 2 * T_explicit[j][i + 1][l]) / (self.dr ** 2)) +
                                                                         ((T_explicit[j][i][l + 1] + T_explicit[j][i][l - 1] - 2 * T_explicit[j][i][l]) / (self.dz ** 2)))
                    else:
                        if l == self.nz_ex - 1:  # r middle and z end
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + \
                                             alpha * self.dt * (((T_explicit[j][i + 1][l] - T_explicit[j][i - 1][l]) / (2 * self.r[i] * self.dr)) +
                                                                ((T_explicit[j][i + 1][l] + T_explicit[j][i - 1][l] - 2 * T_explicit[j][i][l]) / (self.dr ** 2)) +
                                                                ((T_explicit[j][i][l] + T_explicit[j][i][l - 2] - 2 * T_explicit[j][i][l - 1]) / (self.dz ** 2)))
                        elif l == 0:  # r middle and z start
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + (
                                        (2 * self.dt) / (rho * cp * self.dr * self.dz)) * \
                                                      ((k * self.dz * (T_explicit[j][i + 1][l] - T_explicit[j][i][l]) / (2 * self.r[i] * self.dr)) +
                                                       (k * self.dz * (T_explicit[j][i][l] - T_explicit[j][i - 1][l]) / (2 * self.r[i] * self.dr)) +
                                                       (k * self.dr * (T_explicit[j][i][l + 1] - T_explicit[j][i][l]) / self.dz) +
                                                       self.h_m * self.dr * (self.Tc - T_explicit[j][i][0]))
                        else:  # r middle and z middle
                            T_explicit[j + 1][i][l] = T_explicit[j][i][l] + alpha * self.dt * (((T_explicit[j][i + 1][l] - T_explicit[j][i - 1][l]) /
                                                                  (2 * self.r[i] * self.dr)) + ((T_explicit[j][i + 1][l] +
                                                                  T_explicit[j][i - 1][l] - 2 *T_explicit[j][i][l]) / (self.dr ** 2)) +
                                                                  ((T_explicit[j][i][l + 1] + T_explicit[j][i][l - 1] - 2 * T_explicit[j][i][l]) / (self.dz ** 2)))
        return T_explicit

if __name__ == "__main__":
    # CASING ANALYSIS
    material_case = 'aluminium6061t6'
    material_liner = 'epdm'
    insulatior_thk = 0.003
    case_thk = 0.00555
    r_steps = 20
    hm = 1295
    Tc = 1600
    Ta = 297
    lenght = 1.42
    burn_time = 5
    t_steps = 5000
    type = 2
    motor = motor_case(material_case, material_liner, insulatior_thk, case_thk, r_steps, hm, Tc, Ta, lenght, burn_time, t_steps, type)
    print(motor.T)

    # BULKHEAD ANALYSIS
    """
    t_steps_ex = 500
    z_steps_ex = 20
    r_steps_ex = 500
    insulator_thk_ex = 0.003
    z_thk_ex = 0.01
    r_ex_ex = 0.11
    r_in_ex = 0.12
    bulk = bulkhead(material_case, material_liner, insulator_thk_ex, case_thk, r_steps_ex, hm, Tc, Ta, lenght, burn_time,
                    t_steps_ex, z_thk_ex, r_ex_ex, r_in_ex)
    print(bulk.T)
    """