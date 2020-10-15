import numpy as np
from materials import *
from gui import *

class motor_case():
    def __init__(self, material_case, material_liner,  insulator_thk, case_thk, r_steps, hm, Tc, Ta, motor_lenght,  burn_time, t_steps):
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
        self.alpha_case =  float(self.k_case / (self.cp_case * self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))
        self.T = self.run_analysis()

    def run_analysis(self):
        # Radial coordinate
        self.r_1 = 0.11  # Radial position of insulation beginning [m]
        self.r_2 = self.r_1 + self.insulator_thk  # Radial position of insulation/casing interface [m]
        self.r_3 = self.r_2 + self.case_thk  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr  # r variation [m]
        self.nr = self.sectionsr + 1  # Nº of radial points

        # Height coordinate
        self.dz = self.dr  # z variation [m]
        self.sectionsz = int(np.ceil(self.z / self.dz))  # Nº of sections beewtween 0 and z
        print(self.sectionsz)
        self.nz = self.sectionsz + 1  # Nº of z points

        # Time
        self.dt = self.t / self.nt  # Time variation [s]

        # we have to set the initial and boundary contitions
        self.r = [self.r_1 + i * self.dr for i in range(self.nr)]
        self.z = [i * self.dz for i in range(self.nz)]

        A = self.create_Amatrix()

        # Creating the T matrix
        T = np.zeros([self.nt, self.nr])

        # Setting initial values for the T matrix
        for i in range(self.nr):
            T[0][i] = self.Ta

        # Getting the inverse matrix of A
        A_inverse = np.linalg.inv(A)

        # Calculating the T matrix
        for i in range(self.nt - 1):
            T[i][0] = T[i][0] + 2 * self.dt * self.h_m * (self.Tc - T[i][0]) / (self.rho_insulator * self.cp_insulator * self.dr)
            T[i + 1] = (np.matmul(A_inverse, T[i]))

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


if __name__ == "__main__":
    import sys
    material_case = 'aluminium6061t6'
    material_liner = 'epdm'
    insulatior_thk = 0.003
    case_thk = 0.00555
    r_steps = 500
    hm = 1295
    Tc = 1600
    Ta = 297
    lenght = 1.42
    burn_time = 5
    t_steps = 500
    motor = motor_case(material_case, material_liner, insulatior_thk, case_thk, r_steps, hm, Tc, Ta, lenght, burn_time, t_steps)
    print(motor.T)