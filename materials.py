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
nbr = insulation(1000, 0.24, 1350)
paper = insulation(850, 0.18, 1340)
fiberglass = insulation(1900, 0.3, 1250)
phenolic = insulation(1356, 0.29, 1400)


aluminium6061t6 = case(2700, 167, 896)
stainlessteel304 = case(7900, 19, 500)
stainlessteel316 = case(8000, 16.3, 500)
steel1010 = case(7870, 49.8, 650)
carbonfiber = case(1850, 0.8, 1000)




def case_selector(case: str):
    if case == 'Aluminium 6061-T6':
        rho_case, k_case, cp_case = aluminium6061t6.rho_case, aluminium6061t6.k_case, aluminium6061t6.cp_case
    if case == 'Stainless Steel 304':
        rho_case, k_case, cp_case = stainlessteel304.rho_case, stainlessteel304.k_case, stainlessteel304.cp_case
    if case == 'Stainless Steel 316':
        rho_case, k_case, cp_case = stainlessteel316.rho_case, stainlessteel316.k_case, stainlessteel316.cp_case
    if case == 'Steel 1010':
        rho_case, k_case, cp_case = steel1010.rho_case, steel1010.k_case, steel1010.cp_case
    if case == 'Carbon Fiber':
        rho_case, k_case, cp_case = carbonfiber.rho_case, carbonfiber.k_case, carbonfiber.cp_case

    return rho_case, k_case, cp_case

def insulator_selector(liner: str):
    if liner == 'EPDM':
        rho_insulator, k_insulator, cp_insulator = epdm.rho_insulator, epdm.k_insulator, epdm.cp_insulator
    if liner == 'NBR':
        rho_insulator, k_insulator, cp_insulator = nbr.rho_insulator, nbr.k_insulator, nbr.cp_insulator
    if liner == 'Paper':
        rho_insulator, k_insulator, cp_insulator = paper.rho_insulator, paper.k_insulator, paper.cp_insulator
    if liner == 'FiberGlass':
        rho_insulator, k_insulator, cp_insulator = fiberglass.rho_insulator, fiberglass.k_insulator, fiberglass.cp_insulator
    if liner == 'Phenolic Paper':
        rho_insulator, k_insulator, cp_insulator = phenolic.rho_insulator, phenolic.k_insulator, phenolic.cp_insulator

    return rho_insulator, k_insulator, cp_insulator