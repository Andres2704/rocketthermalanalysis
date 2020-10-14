import math
# Liner properties
k_insulator = 0.2                                                                           # Insulation thermal conductivity [W/m-K]
rho_insulator = 860                                                                         # Insulation density [kg/m3]
cp_insulator = 2000                                                                         # Insulation specific heat capacity [J/kg-K]
alpha_insulator = float(k_insulator / (cp_insulator * rho_insulator))                       # Insulation thermal diffusivity [m^2/s]
# ------------------------------------------------
# Case properties
k_case = 167                                                                                # Casing thermal conductivity [W/m-K]
rho_case = 2700                                                                             # Casing density [kg/m3]
cp_case = 896                                                                               # Casing specific heat capacity [J/kg-K]
alpha_case = float(k_case / (cp_case * rho_case))                                           # Casing thermal diffusivity [m^2/s]
# -----------------------------------------------
# Convection Coefficient
h_m = 1295                                                                                  # Average convection coefficient [W/m^2-K]
#-------------------------------------------------
# Temperature
Ta = 297                                                                                    # Ambient temperature in Tatuí-SP-Brazil in the month of August [K]
Tc = 1600                                                                                   # Chamber temperature [K]
# ------------------------------------------------
# Inputs
# Radial coordinate
insulator_thk = 0.003                                                                       # Insulation thickness [m]
case_thk = 0.00555                                                                          # Case thickness |VERIFY THIS VALUEEEE| [m]
r_1 = 0.11                                                                                  # Radial position of insulation beginning [m]
r_2 = r_1 + insulator_thk                                                                   # Radial position of insulation/casing interface [m]
r_3 = r_2 + case_thk                                                                        # Radial position of casing end [m]
sectionsr = 500                                                                             # Nº of sections between r_1 and r_3
dr = (r_3 - r_1) / sectionsr                                                                # r variation [m]
nr = sectionsr + 1                                                                          # Nº of radial points

# Height coordinate
z = 1.42                                                                                    # Height of combustion chamber [m]
dz = dr                                                                                     # z variation [m]
sectionsz = math.ceil(z / dz)                                                               # Nº of sections beewtween 0 and z
nz = sectionsz + 1                                                                          # Nº of z points

# Time
t = 5                                                                                      # Burn time [s]
nt = 500
dt = t / nt                                                                                 # Time variation [s]

# Stability constants
Fo_insulator = alpha_insulator * (dt / (dr ** 2))                                           # Insulation Fourier number
Fo_case = alpha_case * (dt / (dr ** 2))                                                     # Casing Fourier number

Bi_insulator = h_m * (dr) / k_insulator                                                     # Insulation Biot number
Bi_case = h_m * (dr) / k_case                                                               # Casing Biot number