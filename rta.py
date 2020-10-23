from PyQt5 import QtWidgets
from gui import Ui_RocketThermalAnalysis

import sys, math, materials
import numpy as np
import datetime, os
import pandas as pd
import matplotlib.pyplot as plt


class myWindows(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_RocketThermalAnalysis()  # Calling the main windows of the UI
        self.ui.setupUi(self)
        self.case_temperature = []  # Preallocating the case temperature object
        self.bulkhead_temperature = []  # Preallocatin the bulkhead temperature object
        self.run_implicit_button = self.ui.run_implicit.clicked.connect(
            self.run_implicit_function)  # Run implicit method for case
        self.run_explicit_button = self.ui.run_explicit.clicked.connect(
            self.run_explicit_function)  # Run explicit method for case
        self.run_bulk_button = self.ui.run_bulk.clicked.connect(
            self.run_bulkhead_function)  # Run bulkhead explicit method
        self.generate_out_case_button = self.ui.generate_output.clicked.connect(
            self.generate_output_function)  # Generate output file for case analysis
        self.generate_out_bulk_button = self.ui.generate_output_2.clicked.connect(
            self.generate_output_function_bulk)  # Generate output file for bulkhead analysis
        self.generate_graph_case = self.ui.generate_graphs.clicked.connect(
            self.generate_graph_case_function)  # Generate graphs for case analysis
        self.generate_graph_bulk = self.ui.generate_graphs_2.clicked.connect(
            self.generate_graph_bulk_function)  # Generate graphs for bulkhead analysis

    def display_errors(self, text, message):
        """
        This method is for display an error message box to the user, widely used along the code

        :param text: The header text which is shown in the message box

        :param message: The main message

        :return: A message box
        """
        error = QtWidgets.QMessageBox()
        error.setIcon(error.Critical)
        error.setText(text)
        error.setInformativeText(message)
        error.setWindowTitle('Error')
        error.exec_()

    def display_loadmessage(self, text, message):
        """
         This method is for display an information message box to the user, widely used along the code

        :param text: The header text which is shown in the message box

        :param message: The main message

        :return: A message box
        """
        error = QtWidgets.QMessageBox()
        error.setIcon(error.Information)
        error.setText(text)
        error.setInformativeText(message)
        error.setWindowTitle('Atention')
        error.exec_()

    def display_loadwarning(self, text, message):
        """
         This method is for display a warning message box to the user, widely used along the code

        :param text: The header text which is shown in the message box

        :param message: The main message

        :return: A message box
        """
        error = QtWidgets.QMessageBox()
        error.setIcon(error.Warning)
        error.setText(text)
        error.setInformativeText(message)
        error.setWindowTitle('Warning')
        error.exec_()

    def run_implicit_function(self):
        """
        This method is used when the user press to analyse an IMPLICIT numerical method for the rocket motor cases

        :return: Returns nothing else than 1. The main purpose for the function is to declare an object of the
        motor case class and define a property of myWindows classes with the resultant analysis
        """
        try:
            # Saving the data provided by the user
            material_liner = self.ui.material_insulator.text()
            material_case = self.ui.material_case.text()
            insulator_thk = self.ui.insulator_thk.text()
            case_thk = self.ui.case_thk.text()
            ri = self.ui.ri.text()
            burn_time = self.ui.burn_time.text()

            t_steps = self.ui.t_steps.text()
            if int(t_steps) < 10:  # Setting a minimum value to run the analysis
                self.display_errors('Minimum value', 'Please put at least 10 points for time step.')
                return 0

            r_steps = self.ui.r_steps.text()
            if int(r_steps) < 10:  # Setting a minimum value to run the analysis
                self.display_errors('Minimum value', 'Please put at least 10 points for radial step.')
                return 0

            if int(t_steps)>500 or int(r_steps)>1000:
                self.display_loadwarning('Output', 'This time OR radial step will generate a weighty output file. [Ok to run the analysis]')

            hm = self.ui.hm.text()
            Ta = self.ui.Ta.text()
            Tc = self.ui.Tc.text()

            # Verify if some field is empty
            self.data = [material_case, material_liner, insulator_thk, case_thk, ri, t_steps, burn_time, r_steps, hm,
                         Ta, Tc]
            for i in self.data:
                if i == '':
                    self.display_errors('Empty field', 'There is some empty field, please fill it and run again')
                    return 0

            self.method = 'Implicit'  # This will be used when we generate the graphical resources

            # Calling the solution method
            self.case_temperature = motor_case(material_case, material_liner, insulator_thk, case_thk,
                                               ri, t_steps, burn_time, r_steps,
                                               hm, Ta, Tc,
                                               1)  # The last parameter is to define the type of analysis | 1-Implicit | 2 - Explicit
            return 1
        except:
            # Display an error in case the of failure of the analysis
            self.display_errors('Memory error',
                                'As a result from your data a memory error has result, please check it.')
            return 1

    def run_explicit_function(self):
        """
         This method is used when the user press to analyse an EXPLICIT numerical method for the rocket motor cases

         :return: Returns nothing else than 1. The main purpose for the function is to declare an object of the
         motor case class and define a property of myWindows classes with the resultant analysis
         """
        try:
            # Saving the text provided by the user
            material_liner = self.ui.material_insulator_2.text()
            material_case = self.ui.material_case_2.text()
            insulator_thk = self.ui.insulator_thk_2.text()
            case_thk = self.ui.case_thk_2.text()
            ri = self.ui.ri_2.text()
            burn_time = self.ui.burn_timeexp.text()
            hm = self.ui.hm_2.text()
            Ta = self.ui.Ta_2.text()
            Tc = self.ui.Tc_2.text()

            t_steps = self.ui.t_steps_exp.text()
            if int(t_steps) < 10:  # Setting a minimum value to run the analysis
                self.display_errors('Minimum value', 'Please put at least 10 points for time step.')
                return 0

            r_steps = self.ui.r_steps_2.text()
            if int(r_steps) < 10:  # Setting a minimum value to run the analysis
                self.display_errors('Minimum value', 'Please put at least 10 points for radial step.')
                return 0

            # Verify if some field is empty
            self.data = [material_case, material_liner, insulator_thk, case_thk, ri, t_steps, burn_time, r_steps, hm,
                         Ta, Tc]
            for i in self.data:
                if i == '':
                    self.display_errors('Empty field', 'There is some empty field, please fill it and run again')
                    return 0

            self.method = 'Explicit'
            # Calling the solution method
            self.case_temperature = motor_case(material_case, material_liner, insulator_thk, case_thk,
                                               ri, t_steps, burn_time, r_steps,
                                               hm, Ta, Tc, 2)
            return 1
        except:
            self.display_errors('Memory error',
                                'As a result from your data a memory error has result, please check it.')
            return 1

    def run_bulkhead_function(self):
        """
         This method is used when the user press to analyse a 2D EXPLICIT numerical method for the rocket motor bulkhead

         :return: Returns nothing else than 1. The main purpose for the function is to declare an object of the
         motor case class and define a property of myWindows classes with the resultant analysis
         """
        try:
            # Saving the text provided by the user
            material_liner = self.ui.material_insulator_4.text()
            material_case = self.ui.material_bulk.text()
            insulator_thk = self.ui.insulator_thk_4.text()
            bulkhead_thk = self.ui.z.text()
            r_inner = self.ui.ri_4.text()
            r_outer = self.ui.radius_outer.text()
            burn_time = self.ui.burn_time_4.text()
            r_steps = self.ui.r_steps_4.text()
            hm = self.ui.hm_3.text()
            Ta = self.ui.Ta_3.text()
            Tc = self.ui.Tc_3.text()

            r_steps = self.ui.r_steps_4.text()
            if int(r_steps) < 10:  # Setting a minimum value to run the analysis
                self.display_errors('Minimum value', 'Please put at least 10 points for radial step.')
                return 0

            # Verify if some field is empty
            data = [material_case, material_liner, insulator_thk, bulkhead_thk, r_inner, r_outer, burn_time, r_steps,
                    hm, Ta, Tc]
            for i in data:
                if i == '':
                    self.display_errors('Empty field', 'There is some empty field, please fill it and run again')
                    return 0

            # Calling the solution method
            self.bulkhead_temperature = bulkhead(material_case, material_liner, insulator_thk, bulkhead_thk, r_inner,
                                                 r_outer, burn_time, r_steps, hm, Ta, Tc)

        except:
            self.display_errors('Memory error',
                                'As a result from your data a memory error has result, please check it.')
            return 1

    def generate_output_function(self):
        """
         The purpose of this method is to generate an output file with the results of analysis for the case

         :return: Returns nothing else than 1 and an output file in the same directory that this code is running
         """
        # Verifying if some analysis was generated
        if self.case_temperature == []:
            self.display_errors('Error', 'You did not run any analysis, please run it to generate a graph')
        else:
            # Setting up the dataframe where the resultant analysis will be saved
            data = pd.DataFrame(data=self.case_temperature.T)
            if len(data[0]) <= 1000 and self.method=='Implicit':
                pd.set_option('display.max_rows', None)
            elif self.method=='Explicit':
                # Generating the dataframe for t=0%:100% with a step of 10%
                factor = len(data[0])
                fractions = np.array([i/10 for i in range(0,10)])
                indexes = list(set(np.ceil(fractions*factor)))
                indexes= [int(i) for i in indexes]
                indexes.sort()
                new_data = [list(data.iloc[i]) for i in indexes]
                new_data = new_data + [list(data.iloc[factor-1])]
                data = pd.DataFrame(new_data, index=[str(i*10)+'%' for i in range(0,11)])

            pd.set_option('display.width', None)
            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_colwidth', None)

            # Generating the resultant file
            date = datetime.datetime.now()
            information = '''Output file generated by Rocket Thermal Analysis\n\nType of analysis: ''' + 'Case temperature distribution' + '''
            \nAnalysis Information:\n Case material: ''' + str(self.data[0]) + '''\n Insulator Material: ''' + str(
                self.data[1]) + '''\n Insulator thickness: ''' + str(self.data[2]) + ''' m\n Case thickness: ''' + str(
                self.data[3]) + ''' m\n Inner cylinder radius: ''' + \
                          str(self.data[4]) + ''' m\n Time Steps: ''' + str(self.data[5]) + '''\n Burn time: ''' + str(
                self.data[6]) + ''' s\n Radial Steps: ''' + \
                          str(self.data[7]) + '''\n Convection coefficient: ''' + str(
                self.data[8]) + ''' W/m2-K\n Initial case temperature: ''' + str(
                self.data[9]) + ''' K\n Combustion Temperature: ''' + str(self.data[10]) + ''' K
            \nDate:''' + str(date) + '''
            \nMethod:''' + self.method + '''
            \nThe following dataframe is organized by Time step x Radial position \n\n''' + str(data)

            cwd = os.getcwd()  # Getting the current directory
            filename = '\RTAout' + str(date.strftime('%f')) + '.txt'
            directory = cwd + filename
            file = open(directory, 'a', encoding='ISO-8859-1')
            file.write(information)
            file.close()

            self.display_loadmessage('Output file generated',
                                     'The output file has been generated and it was saved with the .exe file')
            return 1

    def generate_output_function_bulk(self):
        self.display_loadmessage('Unavailable', 'This feature is unavailable at this moment')

    def generate_graph_case_function(self):
        """
         This method is used when the user to generate the graphs of the analysis for the case

         :return: Return three different graphs, are they:
                  - Temperature distribution along the insulator and case
                  - Temperature distribution along the case at the final time
                  - Temperature distribution along the case at time 0, 0.25t, 0.5t, 0.75t and t (t = Burn time)
         """
        if self.case_temperature == []:
            self.display_errors('Error', 'You did not run any analysis, please run it to generate a graph')
        else:
            # Finding the case section
            r_a = np.array(self.case_temperature.r)
            locals = [r_a >= self.case_temperature.r_2]
            r_a = r_a[tuple(locals)]

            # Plot along the liner and the case section
            plt.figure(1, figsize=(10, 20))
            plt.plot(self.case_temperature.r, self.case_temperature.Values, 'b--', linewidth=2,
                     label='Insulator section')
            plt.plot(r_a, self.case_temperature.Values[self.case_temperature.r.index(r_a[0])::], 'r-+', linewidth=2,
                     label='Case section')
            plt.title('Temperature distribution - Insulator thickness: ' + str(
                self.case_temperature.insulator_thk * 1000) + ' mm (t = ' + str(self.case_temperature.t) + 's)')
            plt.xlabel('Radial position [m]')
            plt.ylabel('Temperature [K]')
            plt.legend(title='Legend')

            # Plot only case section
            plt.figure(2, figsize=(10, 20))
            plt.plot(r_a, self.case_temperature.Values[self.case_temperature.r.index(r_a[0])::], 'k', linewidth=2,
                     label='Case section', marker="D", markerfacecolor='r', markeredgecolor='r', markersize=4)
            plt.title('Temperature distribution along the case - Insulator thickness: ' + str(
                self.case_temperature.insulator_thk * 1000) + ' mm (t = ' + str(self.case_temperature.t) + 's)')
            plt.xlabel('Radial position [m]')
            plt.ylabel('Temperature [K]')
            plt.legend(title='Legend')

            # Plot multiple time temperature distribution
            metric = self.case_temperature.nt
            T0 = self.case_temperature.T[0, self.case_temperature.r.index(r_a[0])::]  # T at time zero
            T1 = self.case_temperature.T[round(metric / 4), self.case_temperature.r.index(r_a[0])::]  # T at 25% of t
            T2 = self.case_temperature.T[round(metric / 2), self.case_temperature.r.index(r_a[0])::]  # T at 50% of t
            T3 = self.case_temperature.T[round(metric * (3 / 4)),
                 self.case_temperature.r.index(r_a[0])::]  # T at 75% of t
            T4 = self.case_temperature.T[metric - 1, self.case_temperature.r.index(r_a[0])::]  # T at 100% of t

            plt.figure(3, figsize=(10, 20))
            plt.plot(r_a, T0, label='t = ' + str(self.case_temperature.t * 0))
            plt.plot(r_a, T1, label='t = ' + str(self.case_temperature.t * 0.25))
            plt.plot(r_a, T2, label='t = ' + str(self.case_temperature.t * 0.50))
            plt.plot(r_a, T3, label='t = ' + str(self.case_temperature.t * 0.75))
            plt.plot(r_a, T4, label='t = ' + str(self.case_temperature.t * 1))
            plt.title('Temperature distribution along the case - Insulator thickness: ' + str(
                self.case_temperature.insulator_thk * 1000) + ' mm')
            plt.xlabel('Radial position [m]')
            plt.ylabel('Temperature [K]')
            plt.legend(title='Legend')

            plt.show()

    def generate_graph_bulk_function(self):
        """
         This method is used when the user to generate the graphs of the analysis for the case

         :return: Return three different graphs, are they:
                  - Temperature distribution along the insulator and case
                  - Temperature distribution along the case at the final time
                  - Temperature distribution along the case at time 0, 0.25t, 0.5t, 0.75t and t (t = Burn time)
         """

        if self.bulkhead_temperature == []:
            self.display_errors('Error', 'You did not run any analysis, please run it to generate a graph')
        else:
            plt.plot(self.bulkhead_temperature.bulkhead_points, self.bulkhead_temperature.T_final, 'r-+',
                     label='T = ' + str(self.bulkhead_temperature.t))
            plt.title('Temperature distribution - Insulator thickness: ' + str(
                self.bulkhead_temperature.insulator_thk * 1000) + ' mm')
            plt.xlabel('Vertical position of the bulkhead [m]')
            plt.ylabel('Temperature [K]')
            plt.legend(title='Legend')
            plt.show()


class motor_case(myWindows):
    """
    This class is an inherit class from myWindows, this is because we need to some method which are in myWindows class (i.e
    display message error). Within this class are all the calculation of explicit and implicit method for the motor cases.
    """

    def __init__(self, material_case, material_liner, insulator_thk, case_thk, ri, t_steps, burn_time, r_steps, hm, Ta,
                 Tc, type):
        """
        This is a constructor method of motor_case class

        :param material_case: Material of the cases (Format=str)
        :param material_liner: Material of the insulator(Format=str)
        :param insulator_thk: Insulator thickness (Format=float)
        :param case_thk: Case thickness (Format=float)
        :param ri: Inner cylinder radius (Format=float)
        :param t_steps: Time steps (Format=int)
        :param burn_time: Burning time (Format=float)
        :param r_steps: Radial steps (Format=int)
        :param hm: Heat convective coefficient (Format=float)
        :param Ta: Initial temperature of the materials (Format=float)
        :param Tc: Chamber temperature (Format=float)
        :param type: Type of analysis (Format=int)

        """

        super(motor_case, self).__init__()
        self.insulator_thk = float(insulator_thk)
        self.case_thk = float(case_thk)
        self.r_1 = float(ri)
        self.sectionsr = int(r_steps)
        self.h_m = float(hm)
        self.Tc = float(Tc)
        self.Ta = float(Ta)
        self.t = float(burn_time)
        self.nt = int(t_steps)

        # Importing material data and calculating the conductivitty coefficient of each one of them
        self.rho_case, self.k_case, self.cp_case = materials.case_selector(material_case)
        self.rho_insulator, self.k_insulator, self.cp_insulator = materials.insulator_selector(material_liner)
        self.alpha_case = float(self.k_case / (self.cp_case * self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))
        if type == 1:
            self.T, self.Values = self.run_analysis_case()
        elif type == 2:
            self.T, self.Values = self.run_analysis_explicit()

    def run_analysis_case(self):
        """
         This method is used to calculate the temperature distribution on the motor cases by an implicit method

        :return: Temperature distribution in matrix form and the temperature distribution along the outer case side
        """

        # Radial coordinate
        self.r_2 = self.r_1 + self.insulator_thk  # Radial position of insulation/casing interface [m]
        self.r_3 = self.r_2 + self.case_thk  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr  # r variation [m]
        self.nr = self.sectionsr + 1  # Nº of radial points

        # Time
        self.dt = self.t / self.nt  # Time variation [s]

        # we have to set the initial and boundary conditions
        self.r = [self.r_1 + i * self.dr for i in range(self.nr)]  # Defining the r vector

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
            T_implicit[i][0] = T_implicit[i][0] + (self.h_m * 2 * self.dt * (self.Tc - T_implicit[i][0])) / (
                        self.rho_insulator * self.cp_insulator * self.dr)
            T_implicit[i + 1] = (np.matmul(A_inverse, T_implicit[i]))

        values = T_implicit[self.nt - 1, ::]
        self.display_loadmessage('Analysis finished',
                                 'The analysis has been finished, now you can generate your output file and graphical resources')
        return T_implicit, values

    def run_analysis_explicit(self):
        """
         This method is used to calculate the temperature distribution on the motor cases by an explicit method

        :return: Temperature distribution in matrix form and the temperature distribution along the outer case side
        """

        # Radial coordinate
        self.r_2 = self.r_1 + self.insulator_thk  # Radial position of insulation/casing interface [m]
        self.r_3 = self.r_2 + self.case_thk  # Radial position of casing end [m]
        self.dr = (self.r_3 - self.r_1) / self.sectionsr  # r variation [m]
        self.nr = self.sectionsr + 1  # Nº of radial points

        # Time
        self.dt = self.t / self.nt

        # we have to set the initial and boundary conditions
        self.r = [self.r_1 + i * self.dr for i in range(self.nr)]  # Defining the r vector

        d_insulator = self.alpha_insulator * (self.dt / (self.dr ** 2))  # Insulation Fourier number
        d_case = self.alpha_case * (self.dt / (self.dr ** 2))  # Casing Fourier number

        T = np.zeros([self.nt, self.nr])
        for i in range(self.nr):
            T[0][i] = self.Ta

        # Verifying if the method will converge
        if d_insulator > 0.5 or d_case > 0.5:
            message = 'Stability criteria was not achieved, the following values must be less than 0.5' + '\nFourier number - case: ' + str(
                d_case) + \
                      '\nFourier number insulator: ' + str(d_insulator) + \
                      '\n\nTry to decrease the radial steps or increase time steps'
            self.display_errors('Convergence error', message)

            return 0, 0
        else:
            self.display_loadmessage('Running analysis...',
                                     'The analysis will be starting now, it can be time consuming, wait it. [Press ok to start]')

        # Calculating the T-matrix
        for j in range(self.nt - 1):
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
        values = T[self.nt - 1, ::]
        self.display_loadmessage('Analysis finished',
                                 'The analysis has been finished, now you can generate your output file and graphical resources')
        return T, values

    def create_Amatrix(self):
        """
         This method is used to the tridiagonal matrix to solve by implicit method

        :return: Tridiagonal matrix
        """

        # Creating the A matrix
        A = np.zeros([self.nr, self.nr])

        # Building the A matrix
        A[0][0] = 1 + (2 * self.dt * self.k_insulator) / (self.rho_insulator * self.cp_insulator * self.dr ** 2)
        A[0][1] = -self.k_insulator * 2 * self.dt / (self.dr * self.dr * self.rho_insulator * self.cp_insulator)
        A[self.nr - 1][self.nr - 1] = 1 - self.alpha_case * self.dt / (
                    self.r[self.nr - 1] * self.dr) - self.alpha_case * self.dt / (self.dr ** 2)
        A[self.nr - 1][self.nr - 2] = self.alpha_case * self.dt / (
                    self.r[self.nr - 1] * self.dr) + 2 * self.alpha_case * self.dt / (self.dr ** 2)
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


class bulkhead(myWindows):
    """
    This class is an inherit class from myWindows, this is because we need to some method which are in myWindows class (i.e
    display message error). Within this class are all the calculation of explicit and implicit method for the motor cases.
    """

    def __init__(self, material_case, material_liner, insulator_thk, bulkhead_thk, radius_inner, radius_outer,
                 burn_time, r_steps, hm, Ta, Tc):
        """
        This is a constructor method of bulkhead class

        :param material_case: Material of the cases (Format=str)
        :param material_liner: Material of the insulator(Format=str)
        :param insulator_thk: Insulator thickness (Format=float)
        :param bulkhead_thk: Bulkhead thickness (Format=float)
        :param radius_inner: Inner radius of the bulkhead within the case wall(Format=float)
        :param radius_inner: Outer radius of the bulkhead no considering the case (Format=float)
        :param burn_time: Burning time (Format=float)
        :param r_steps: Radial steps (Format=int)
        :param hm: Heat convective coefficient (Format=float)
        :param Ta: Initial temperature of the materials (Format=float)
        :param Tc: Chamber temperature (Format=float)

        """
        super(bulkhead, self).__init__()

        # Explicit attributes for bulkhead
        self.rho_case, self.k_case, self.cp_case = materials.case_selector(material_case)  # Case properties
        self.rho_insulator, self.k_insulator, self.cp_insulator = materials.insulator_selector(
            material_liner)  # Insulator properties
        self.h_m = float(hm)
        self.insulator_thk = float(insulator_thk)
        self.z = float(bulkhead_thk)
        self.radius_inner = float(radius_inner)
        self.radius_outer = float(radius_outer)
        self.t = float(burn_time)
        self.sectionsr = int(r_steps)
        self.Ta = float(Ta)
        self.Tc = float(Tc)
        self.alpha_case = float(self.k_case / (self.cp_case * self.rho_case))
        self.alpha_insulator = float(self.k_insulator / (self.cp_insulator * self.rho_insulator))

        self.T_bulkhead, self.T_final = self.run_analysis_bulkhead()

    def run_analysis_bulkhead(self):
        """
         This method is used to calculate the 2D temperature distribution on the motor bulkhead by an explicit method

        :return: Temperature distribution in matrix form and the temperature distribution along vertical direction of
        the bulkhead
        """

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
        z = [z_1 + i * dz for i in range(nz)]

        T = np.zeros([nt, nr, nz])
        for i in range(nr):
            for j in range(nz):
                T[0][i][j] = self.Ta

        self.display_loadmessage('Running analysis...',
                                 'The analysis will be starting now, it can be time consuming since is an explicit bidimensional method, please be patient. \n[Press ok to start]')
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
                            T[i + 1][j][l] = (2 * self.h_m * dt * self.Tc) / (rho * cp * dz) + (
                                    1 + (k * dt) / (rho * cp * (dr) ** 2) - (k * dt) / (rho * cp * (dz) ** 2) - (
                                    2 * self.h_m * dt) / (rho * cp * dz)) * T[i][j][l] + (
                                                     (-k * dt) / (rho * cp * (dr) ** 2)) * T[i][j - 1][l] + (
                                                     (k * dt) / (rho * cp * (dz) ** 2)) * T[i][j][l + 1]
                        elif j == 0:  # z start and r start
                            T[i + 1][j][l] = (2 * self.h_m * dt * self.Tc) / (rho * cp * dz) + (
                                    1 - (k * dt) / (rho * cp * (dr) ** 2) - (k * dt) / (rho * cp * (dz) ** 2) - (
                                    2 * self.h_m * dt) / (rho * cp * dz)) * T[i][j][l] + (
                                                     (k * dt) / (rho * cp * (dr) ** 2)) * T[i][j + 1][l] + (
                                                     (k * dt) / (rho * cp * (dz) ** 2)) * T[i][j][l + 1]
                        else:  # z start and r middle
                            T[i + 1][j][l] = (2 * self.h_m * dt * self.Tc) / (rho * cp * dz) + (
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
        self.display_loadmessage('Analysis finished',
                                 'The analysis has been finished, now you can generate your output file and graphical resources')
        T_final = T[nt - 1, nr - 1, ::]
        self.bulkhead_points = z

        print(T)
        print(T_final)

        return T, T_final


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    application = myWindows()
    application.show()
    sys.exit(app.exec())
