#" IN GOD WE TRUST, ALL OTHERS MUST BRING DATA"
#                                               -W. Edwards Deming
#------------------------------------------------------------------------------
# Copyright 2023 The Gamlab Authors. All Rights Reserved.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#------------------------------------------------------------------------------
''' 
The Scientific experimental simulation library 
-------------------------------------------------------------------------------
Graphen & Advanced Material Laboratory 

it aimes to provide new scientist to use data,simlation, prepared data 
and Artificial intelligence models.

See http://gamlab.aut.ac.ir for complete documentation.
'''
__doc__='''

@author: Ali Pilehvar Meibody (Alipilehvar1999@gmail.com)

                                         888                    888
 .d8888b    .d88b.     88888b.d88b.      888         .d88b.     888
d88P"      d88""88b    888 "888 "88b     888        d88""88b    88888PP
888  8888  888  888    888  888  888     888        888  888    888  888
Y88b.  88  Y88..88PP.  888  888  888     888......  Y88..88PP.  888  888
 "Y8888P8   "Y88P8888  888  888  888     888888888   "Y88P8888  88888888  


@Director of Gamlab: Professor M. Naderi (Mnaderi@aut.ac.ir)    

@Graphene Advanced Material Laboratory: https://www.GamLab.Aut.ac.ir


@Co-authors: 
'''




#Constant-------------------------------------------
pi=3.1415926536
e=2.7182818284
c=299792458
N_a=6.022*10**23
R_J=8.3143
R_Cal=1.987
R_LA=0.08205
Pi=3.14
g=9.81
Tau=Pi*2
e=7.1828
pi = 3.14159
e = 2.71828
g = 9.81                                          # Gravitational constant
h =  6.634 * 10**-34                              # constant pelank
k , b = 1.38 * 10 ** -23                          # J/mjck constant of boltzman       
tau = 6.28
R =0.08206           #unit = ((L*atm) / (mol*k)) # Global constant of gases
pi=3.14
ta=6.28
e=2.718
zeta=1.202
#Golden ratio
phi=1.618
N_a = 6.022*10**23              #molecules/mol (Avogadro's Number)
k = 1.381*10**(-23)             #J/K (Boltzmann's Constant)
alpha = 5.67*10**(-8)           #W/m^2/K^4 (Stefan-Boltzmann Constant)
R = 8.314                       #J/mol/K (Gas costant)
Visible = range(380, 750, 10)   #Visible Wavelength
h = 6.626*10**(-34)             #J.s (Planck's Constant)
c = 2.998*10**8                 #m/s (speed of light)
pi=3.14
R=8.314
Na=6.022*(10**23)
e=2.7182818284590452353602874713527
Boltzmann_Constant=1.380649*(10**(-23))
Conductivity_PPY=105
Conductivity_pT=33.7
Conductivity_P3HT=2.4
Conductivity_PLN=10

Kw = 1*10**-14  # Autoionization constant for water at 25°C
Fe_Tm_Alpha=910
'دمای ذوب آهن در فاز آلفا است'

Fe_Tm_Gama=1495
'دمای ذوب آهن در فاز گاما است'

Fe_Tm_Delta=1539
'دمای ذوب آهن در فاز دلتا است'

Fe_Density=7.87
'چگالی آهن و فولاد'

max_C_inSteel=2.11
'بیشترین درصد کربن در فولاد'
'در واقع بیشتر از این مقدار کربن وارد فولاد شود به چدن تبدیل میشود'

Eutectic_T=1148
'دمای واکنش یوتکتیک در دیاگرام آهن کربن'

Eutectic_Percent=4.3
'درصد کربن در فولاد که در آن واکنش یوتکتیک رخ میدهد'

Eutectoid_T=727
'دمای واکنش یوتکتوئید'

Eutectoid_persent=0.76
'درصد کربن در واکنش یوتکتوئید'


Pi=3.14159265358979323846264338327950288419716939937510
C=299792458
#m/s; speed of light in vacuum
R= 8.314510
#J/(mol*K); molar gas constant
g0=9.81 #m/s2

density_of_Fe=7.87
'Here, density of Fe is reported in g per cubic cm'

density_of_Al=2.7
'Here, density of Al is reported in g per cubic cm'

density_of_Cu=8.96
'Here, density of Cu is reported in g per cubic cm'

melting_point_of_Fe=1538
'Here, melting point of Fe is reported in kelvin'

melting_point_of_Al=660
'Here, melting point of Al is reported in kelvin'

melting_point_of_Cu=1085
'Here, melting point of Cu is reported in kelvin'

thermal_conductivity_coefficient_of_Fe=221
'Here, thermal conductivity coefficient of Fe is reported in W/mK'

thermal_conductivity_coefficient_of_Al=237
'Here, thermal conductivity coefficient of Al is reported in W/mK'

thermal_conductivity_coefficient_of_Cu=385
'Here, thermal conductivity coefficient of Cu is reported in W/mK'



R = 8.314                  #Ideal Gas constant (J/mol.K)
h = 6.626 * (10 ** (-34))  #Plank's constant (J.s)
NA = 6.022 * (10 ** (23))  #Avagadro's constant (molecules per mole)
k = 1.381 * (10 ** (-23))  #Boltzman constant (J/K)
c = 3.00 * (10 ** (8))     #Speed of Light is a Vacuum (m/s)
G = 6.674 * (10 ** (-11))  #Universal Gravitational constant ((m^3)/(kg.(s^2)))
F = 96485                  #Faraday's constant (C/mol)


# Define molecular weights of composite materials
G_Mol_Si = 28.09   #silicium
G_Mol_O = 16.00    #oxygen
G_Mol_Ba = 137.33  #barium
G_Mol_Zr = 91.22   #zirconium
G_Mol_Ti = 47.87   #titanium

# Define lists of composite materials
Resin_Materials = ['TEGDMA', 'BIS_GMA', 'UDMA']
Filler_Materials = ['Silica dioxide', 'Barium silicate', 'Zirconium dioxide']
Colors_Materials = ['Minerals: Titanium oxide', 'Organic: Photoactive pigments']




#variable
k= 534.6
pi= 3.14
F=96485
R=8.314
g=23

F=100
Zf=0.25
Xd=0.95
Xw=0.05
R=5
alpha=3
q=0.7



#[Haniyeh] رشته مهندسی کامپیوتر
#dae in ja ma yekseri aadad sabet taeerif kardim
Avogadro_Number=6.02*10**23
#(mol) adad avogadro

Speed_Of_Light=2.99*10**8 
#(m/s) sorat noor--->c

Electron_Charge=1.6*10**19
#(C) bar electeron

Faraday_Constant=9.64*10**4
#C/mol sabet faraday

Earth_Accel=9.8
#(m/s**2) geranesh zamin

        
Planck_constant=6.626*10**(-34) #js#
R=8.314 #J/mol.K#
Faraday_constant= 96485 #C/mol#
Boltzman_constany= 1.381*10**(-23) #J/K#
Avogardo_Num=6.022*10**(23) #mol^-1#





####constants
s=1###ehtemal tamam pishamadha
a=9.8###shetab dar harekat soghoot azad
Income_Tax_Perecent=0.09
Iran_Stock_Market_Volatility=0.05
Santiago_Bernabeu_Stadium_Crowd=78000

#Gravitational_Acceleration_Constant
g=9.81#(m/s^2)

#Standard_Atmospheric_Pressure
P_0=101325#(Pa)

#Plank's_Constant
h=6.62607015e-34 #(J.s)


#Boltzmann_constant
K=1.380649e-23 #(J/K)

#Avogadro's_constant
N_A=6.022e23 #(without unit) 




#constants
Latent_heat = 1.16 * (10 ** 9)
Gama_sl = 0.132
T_m = 1064 #in centigrade
A = 0.413 #in nano meter
Pi = 3.14
K = 1.38 * (10 ** (-23)) #Boltzmann constant
N_t = 3 * (10 ** 33)


N_a = 6.022*10**23              #molecules/mol (Avogadro's Number)
k = 1.381*10**(-23)             #J/K (Boltzmann's Constant)
alpha = 5.67*10**(-8)           #W/m^2/K^4 (Stefan-Boltzmann Constant)
R = 8.314                       #J/mol/K (Gas costant)
Visible = range(380, 750, 10)   #Visible Wavelength
h = 6.626*10**(-34)             #J.s (Planck's Constant)
c = 2.998*10**8                 #m/s (speed of light)





#Part1 ------------------------------------------------------------------------
#IE_Constants

k=1.96              #k=Z(0.05) Normal Standard Distribution
pf=0.6209           #f(P/F , i% , n) = 0.6209 for i=10% and n=5
ag=3.1936           #f(A/G , i% , n) = 3.1936 for i=18% and n=10
t_stu=0.9277        #T-Student Distribution for x=1.5 and n=30 (degree freedom 29)
d=1.128             #Constant coefficients for n=2



###----------Part1---------------
R=8.314 ##gas constant in J/mol.K
D=2.3*(10**(-5)) ##oxygen diffusion coefficient in water μm2/s
S=28.34 ##standard enthropy for solid alumminuim in J/mol.K
k=1.38*(10**(-23))  ##boltzmann constant for gas particle kinetic energy
Tm=14.025 ## hydrojen melting point in K

#part 1:
vacuum_permeability=1
π=3.14
G=6.674 * 10^-11 
e=2.71828
k=9*10**9

#Building importance factor
I1=1.4
I2=1.2
I3=1
I4=0.8
#Length to width ratio of works
gold=618.1

#first________________________________________________________________________________________________
R_J=8.314 #J/mol.K : The molar gas constant or ideal gas constant.
R_Lit=0.08 #Lit.atm/mol.K
R_cal=1.98 #Cal/mol.K
R_cm3=82.06 #cm^3.atm/mol.K
#second________________________________________________________________________________________________
C_V_J_1=12.471 #J/mol.K :The Monoatomic ideal gas constant-VOLUME specific heat.:: C_V=1.5R
C_V_J_2=20.785 #J/mol.K :The Diatomic ideal gas constant-VOLUME specific heat.:: C_V=2.5R
C_V_J_more=29.099 #J/mol.K :The Polyatomic ideal gas constant-VOLUME specific heat.:: C_V=3.5R
#_________________________________________________________________________________________________
C_V_Cal_1=2.97 #Cal/mol.K :The Monoatomic ideal gas constant-VOLUME specific heat.:: C_V=1.5R
C_V_Cal_2=4.95 #Cal/mol.K :The Diatomic ideal gas constant-VOLUME specific heat.:: C_V=2.5R
C_V_Cal_more=6.93 #Cal/mol.K :The Polyatomic ideal gas constant-VOLUME specific heat.:: C_V=3.5R
#third________________________________________________________________________________________________
C_P_J_1=20.785 #J/mol.K :The Monoatomic ideal gas constant-PRESSURE specific heat.:: C_P=2.5R
C_P_J_2=29.099 #J/mol.K :The Diatomic ideal gas constant-PRESSURE specific heat.:: C_P=3.5R
C_P_J_more=37.413 #J/mol.K :The Polyatomic ideal gas constant-PRESSURE specific heat. :: C_V=4.5R
#_________________________________________________________________________________________________
C_P_Cal_1=4.95 #Cal/mol.K :The Monoatomic ideal gas constant-PRESSURE specific heat.:: C_P=2.5R
C_P_Cal_2=6.93 #Cal/mol.K :The Diatomic ideal gas constant-PRESSURE specific heat.:: C_P=3.5R
C_P_Cal_more=8.91 #Cal/mol.K :The Polyatomic ideal gas constant-PRESSURE specific heat. :: C_V=4.5R
#ّfourth________________________________________________________________________________________________
Pi=3.14 #The number π is a mathematical constant that is the ratio of a circle's circumference to its diameter.
#ّfifth________________________________________________________________________________________________
h=6.62*10**(-34) #kg.m^2/s


R = 8.314                        # Gas constant        --> j/mol*k

a = 6.0232 * (10**23)            # Avogadro constant   --> 1/mol

k = 1.38054 * (10**23)           # Boltzmann constant  --> joules/degree

F = 2306                         # Faraday constant    --> calories/volt*mol

h = 6626068 * (10**34)           # Planck constant     --> J*s

c = 4.18                         # Specific heat capacity of liquid water 

m = 9.1093837015 * (10**(-34))   # Electron mass --> kg


k=1.38054*10**-23 
# k is Boltzman constant and in terms jouls/degree

R=8.314 
# R is gasses constant in terms jouls/degree,mole

F=23060
# F is Faraday constant in terms calories/Volt.mole

a=6.0232*10**23
# a is Avogadro constant in terms /gr.mole

h=6.62607015*10**-34
# h is Planck's constant in terms kg.m2.s-1

G=6.674*10*-11
# G is universal gravitational constant 


Boltzmann_Constant=1.380649*(10**(-23)) #joule per kelvin (J/K)
Avogadro_Constant=6.02214*(10**23) #per moles (mol-1)
Faraday_Constant=96485.3399 #coulombs per mole of electrons (C/mol)
Planck_Constant=6.62607015*(10**(-34)) #joule second (J.s)
Elementary_Charge=1.602176634*(10**(-19)) #coulombs (C)
Light_Speed_Constant=299792458 #meters per second (m/s)


