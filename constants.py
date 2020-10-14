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

