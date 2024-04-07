# Rocket Thermal Analysis
### Code created by:
* Andres Benoit (Aerospace Engineering Undegraduate Student - UFSM) - andres.benoit7@gmail.com
* João Vitor Rohr (Aerospace Engineering Undegraduate Student - UFSM) - rohr_joaovitor@hotmail.com

### Overview
Rocket Thermal Analysis(RTA) is an open-source thermal analysis software developed by Tau Rocket Team for the Latin American Space Challenge Cooperation 2021. The software estimates the temperature distribution along the case and bulkhead of a rocket motor.

*Current Features*:
- Temperature distribution along the case by implicit and explicit methods.
- Temperature distribution along the bulkhead by 2D explicit method.
- Generate multiple graphs and output file with the results.
- Convective heat transfer coefficient calculator.
- Multiple and custom material selector
- Open and save study case with .RTA file (feature added in 07/04/2024)

*Future features*:
- Create a pop-up window to build specific graphs
- .xlsx file with the results
- Coasting phase 

### How to use RTA(Rocket Thermal Analysis):
* You can download the executable file [here](https://mega.nz/file/xJF0xByB#TC003mz6IPvXfgTSPowNxX6tO_Ud2U0WDiWHZQrwk3A).
Once you have opened RTA the main window will be something like this:

![alt text](https://github.com/Andres2704/rocketthermalanalysis/blob/master/images/main.PNG)

The main window is destinated to analyse the 1D temperature distribution along the motor case and for that you have two options of numerical methods. One is the implicit (recommended) and the other is the explicit method, the second one is more time consuming than the first. 

On the second tab you will be able to simulate a 2D temperature distribution along the motor's bulkhead. This is a slow method since it is explicit and bidimensional, so be patient.

All the three analysis can be exported as output txt file which you can easily import to excel or other sheet editor and also you can generate plots of the results. 

For all the mathematical model of the software you can find [here](https://github.com/Andres2704/rocketthermalanalysis/blob/master/Thermal_Analysis_software.pdf).

### Source:

All the code were implemente in Python and the main UI was developed with PyQt5 library and the following libraries were also used:
* Numpy
* Matplotlib
* Pandas
