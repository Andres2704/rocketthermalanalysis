class case():
    def __init__(self, rho, k, cp):
        self.rho_case = rho
        self.k_case = k
        self.cp_case = cp

class insulation():
    def __init__(self, rho, k, cp):
        self.rho_insulator = rho
        self.k_insulator = k
        self.cp_insulator = cp

epdm = insulation(860, 0.2, 2000)
aluminium6061t6 = case(2700, 167, 896)

def case_selector(case: str):
    if case == 'Aluminum 6061 - T6':
        rho_case, k_case, cp_case = aluminium6061t6.rho_case, aluminium6061t6.k_case, aluminium6061t6.cp_case

    return rho_case, k_case, cp_case

def insulator_selector(liner: str):
    if liner == 'EPDM':
        rho_insulator, k_insulator, cp_insulator = epdm.rho_insulator, epdm.k_insulator, epdm.cp_insulator

    return rho_insulator, k_insulator, cp_insulator