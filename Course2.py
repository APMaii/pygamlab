#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 14:42:26 2024

@author: apm


Babak yadegari 1 function
alireza peymani 1 function


"""

'''

#========================
#========================
#========================
#========================
#========================
#========================

#========================
#A1


maryam haj heydari
mojtaba bakhshi
shahin esfandiari
masooome
narges tehrani
atena jalilzadeh
erfan jalilzadeh
narges bagherie
zahra afshar
niloofar kazemkhanloo
nazanin nasiri

shahrokh akbari
hamed taghipoor
parinaz hajiali
hanie salehi
zahra moradi
dental_material
mohammadmahdi yektafallah

pavane seyedpour
zahra moradi
nasim heidarian


'''

class Electric_Charge_Convertor:
    def Convert_Coulomb_to_Electron_volt(self , coulomb):
        electron_volt = coulomb * 6.24e18 
        return electron_volt
    
    def Convert_Electron_volt_to_Coulomb(self , electron_volt):
        coulomb = electron_volt / 6.24e18
        return coulomb
    
coulomb_value = float(input("Enter the value in coulomb: "))
electron_volt_value = float(input("Enter the value in electron volt: "))
convertor = Electric_Charge_Convertor()

converted_electron_volt = convertor.Convert_Coulomb_to_Electron_volt(coulomb_value)
print(f"{coulomb_value} coulomb is equal to {converted_electron_volt} electron volts.")

converted_coulomb = convertor.Convert_Electron_volt_to_Coulomb(electron_volt_value)
print(f"{electron_volt_value}  electron volts value is equal to {converted_coulomb} coulomb")

        

def Convert_Gas_Constant(gas_constant , from_unit , to_unit):
    
    conversion_factors = {
        "J/mol.K": {
            "J/mol.K": 1 ,
            "cal/mol.K": 0.239 , 
            "atm.L/mol.K": 8.314 , 
            "cm^3.atm/mol.K": 82.057 , 
        } , 
        "cal/mol.K": {
            "J/mol.K": 4.184 ,
            "cal/mol.K": 1 , 
            "atm.L/mol.K": 24.205 ,
            "cm^3/mol.K": 239.006 , 
        } ,
        "atm.L/mol.K": {
            "J/mol.K": 0.0821 ,
            "cal/mol.K": 0.042 , 
            "atm.L/mol.K": 1 ,
            "cm^3/mol.K": 10 
        } , 
        "cm^3/mol.K": {
            "J/mol.K": 0.001987 ,
            "cal/mol.K": 0.0042 , 
            "atm.L/mol.K": 0.1 ,
            "cm^3/mol.K": 1
        }
    }

    if from_unit in conversion_factors and to_unit in conversion_factors[from_unit]:
        conversion_factors = conversion_factors[from_unit][to_unit]
        converted_value = gas_constant * conversion_factors
        
        return converted_value
    
    else:
        return "Invalid"
    
gas_constant = float(input("Gas constant: "))
from_unit = input("current unit (J/mol.K , cal/mol.K , atm.L/mol.K , cm^3/mol.K): ")
to_unit = input("the unit to convert (J/mol.K , cal/mol.K , atm.L/mol.K , cm^3/mol.K): ")

result = Convert_Gas_Constant(gas_constant , from_unit , to_unit)

if result != "Invalid units":
    print("The converted gas unit is: {:.4f} {}".format(result , to_unit))
    
else:
    print("Try again.")
    
class Hydrogen_Energy_Level_Calculater:   #defined the class 
    def Calculate_Energy_Level(sef , n):  #This method is defined with input in n to calculate the energy level
        energy = -13.6 / (n ** 2)         #The energy level formula
        return energy                     #Energy value is returned as the output of the function
    
    def validate_input(self , n):         #This method gets input n  wich represent th hydrogen energy number and check if it's valide or not
        if not isinstance(n , int) or n <= 0:
            raise ValueError("Invalid value!")
            
try:
    n = int(input("Enter the level number: "))
    
    calculater = Hydrogen_Energy_Level_Calculater()
    calculater.validate_input(n)
    energy_level = calculater.Calculate_Energy_Level(n)
    
    print(f"THe energy level of the hydrogen with n = {n} is {energy_level} ev.")
    
except ValueError as e:
    print(e)
    
    
    



class Grain_Growth_Calculater:
    def __init__(self):
        
        self.k = self.validate_input("Enter the growt rate coefficient (k): ")
        self.grain_size = self.validate_input("Enter the grain size: ")
        self.time_interval = self.validate_input("Enter the time interval (in second): ")
        
    def validate_input(self , prompt):
        while True:
            user_input = input(prompt)
            try:
                value = float(user_input)
                    
                if value <= 0:
                    print("Value must be +. Try again.")
                        
                else:
                    return value
                    
            except ValueError:
                print("Invalid input.")
                    
    def Growth_Rate_Calculation(self):
        return self.k * self.grain_size / self.time_interval
        
def main():
    grain_growth_calculater = Grain_Growth_Calculater()
    growth_rate = grain_growth_calculater.Growth_Rate_Calculation()
    
    print(f"Grain growth rate is {growth_rate}")
                    
if __name__ == "__main__":
    main()
        
                    
                    
                    
class Kinetic_Energy:
    
    def __init__(self):
        self.mass = None
        self.velocity = None
        
    def User_Inputs(self):
        while True:
            try:
                self.mass = float(input("Enter the mass of the object(kg): "))    #Get the object mass from user (int or float)
                self.velocity = float(input("Enter the velocity of the object(m/s): "))    #Get the velocity of the object (int or float)
                break
            
            except ValueError:
                print("Invalid!")
                
    def Kinetic_Energy_Calculater(self):
        
        kinetic_energy = 0.5 * self.mass * self.velocity ** 2 
        return kinetic_energy
    
    def Result(self , result):
        print("Kinetic Energy is: {:.2f} J".format(result))
        
calculater = Kinetic_Energy()
calculater.User_Inputs()
result = calculater.Kinetic_Energy_Calculater()
calculater.Result(result)
        


import math



radius = float(input("Enter the atomic radius in Nanometer: "))
crystal_structure = input("Enter the crystal structure (SC , BCC , FCC): ")

APF = Atomic_Packing_Factor(radius, crystal_structure)
if APF is not None:
    print(f"The Atomic Packing Factor is {APF}")

R = 8.314                  #Ideal Gas constant (J/mol.K)
h = 6.626 * (10 ** (-34))  #Plank's constant (J.s)
NA = 6.022 * (10 ** (23))  #Avagadro's constant (molecules per mole)
k = 1.381 * (10 ** (-23))  #Boltzman constant (J/K)
c = 3.00 * (10 ** (8))     #Speed of Light is a Vacuum (m/s)
G = 6.674 * (10 ** (-11))  #Universal Gravitational constant ((m^3)/(kg.(s^2)))
F = 96485                  #Faraday's constant (C/mol)




def Carnot_Efficiency(T_c , T_h): 
    
    '''
    Parameters
    ----------
    
    T_hot (float): Hot source temperature (Kelvin)
    T_cold (float): Cold source temperature (Kelvin)
    
    Returns:
    float : Efficiency of the cycle 
    
   '''
    try:
        if T_h <= 0 or T_c <=0:
            raise ValueError("Temperature must be +!")
        efficiency = 1 - (T_c / T_h)
        return efficiency 
        
    except ValueError as e:
        print(f"ERROR: {e}")
        
T_cold = float(input("Enter the cold source temperature: "))
T_hot = float(input("Enter the hot source temperature: "))

efficiency = Carnot_Efficiency(T_cold, T_hot)

if efficiency is not None:
    print(f"Carnot efficiency is {efficiency: .2f}")
        
        
        
        





#1st part
import math
import numpy as np
F=100
Zf=0.25
Xd=0.95
Xw=0.05
R=5
alpha=3
q=0.7
def Mc_Cabe(F,Zf,Xd,Xw,R,alpha,q):
    '''This function is used for Mc-Cabe calculation in Distillation Towers
    F : Feed Rate
    Zf : volatile composition in Feed
    Xd : volatile composition in Distillate
    Xw : volatile composition in Waste
    R : Reflux
    alpha : Volatility coefficient
    q : the quantity of liquid in feed
    Returns the number of tray'''
    W=F*(Xd-Zf)/(Xd-Xw)
    D=F-W
    L_first=R*D
    G_first=(R+1)*D
    L_second=L_first+q*F
    G_second=G_first-(1-q)*F  
    # Rmin Calculation
    guess=0
    n=0
    while n<10000:
        if q==0:
            Left=np.array([Zf/(1-Zf)])
            Right=np.array([alpha*((Xd*(q-1)+Zf*(guess+1))/(((1-Xd)*(q-1))+(1-Zf)*(guess+1)))])
        else:
            Left=np.array([((Xd*q+Zf*guess)/((1-Xd)*q+(1-Zf)*guess))])
            Right=np.array([alpha*((Xd*(q-1)+Zf*(guess+1))/(((1-Xd)*(q-1))+(1-Zf)*(guess+1)))])
        if Left-Right<0:
            Rmin=guess
        else:
            guess=guess+0.001            
        n=n+1
        
    # Nmin Calculation
    Nmin=np.ceil(math.log((Xd*(1-Xw))/(Xw*(1-Xd)),10)/(math.log(alpha,10))-1)
    
    # N Calculation
    # 1st Operating Line
    x=np.array([])
    for i in range(0,101,1):
        a=np.array([i/100])
        x=np.concatenate((x,a))
        
    yeq=np.array([])
    for i in range(0,101,1):
        b=np.array([(alpha*(i/100))/((alpha-1)*(i/100)+1)])
        yeq=np.concatenate((yeq,b))
        
    y_first=np.array([])
    for i in range(0,101,1):
        c=np.array([(((R*(i/100))/(R+1))+((Xd)/(R+1)))])
        y_first=np.concatenate((y_first,c))
        
    y_second=np.array([])
    for i in range(0,101,1):
        d=np.array([L_second*(i/100)/G_second-W*Xw/G_second])
        y_second=np.concatenate((y_second,d))
        
    xfeed=np.array([])
    for i in range(0,101,1):
        if q==1:
            e=np.array([Zf])
            xfeed=np.concatenate((xfeed,e))
        else:
            e=np.array([i/100])
            xfeed=np.concatenate((xfeed,e))
    
    yfeed=np.array([])
    for i in range(0,101,1):
        if q==1:
            m=np.array([i/100])
            yfeed=np.concatenate((yfeed,m))
        else:
            f=np.array([((q*(i/100))/(q-1))-(Zf/(q-1))])
            yfeed=np.concatenate((yfeed,f))
            
        if q==1:
         xcrosspoint=Zf
         ycrosspoint=(L_first*xcrosspoint)/(G_first)+Xd/(R+1)
        else:       
         xcrosspoint=((Zf/(q-1))+(Xd/(R+1)))/((q/(q-1))-(L_first/G_first))
         ycrosspoint=(q*xcrosspoint)/(q-1)-(Zf/(q-1))
    
    Xopt=[]
    x=Xd
    for i in range(0,100):
        if ycrosspoint<x:
            xopt=x/(alpha+x-alpha*x)
            x=(L_first*xopt/G_first)+(Xd/(R+1))
            Xopt=Xopt+[xopt]
            Nfeed=len(Xopt)
        elif x>Xw:
            xopt=x/(alpha+x-alpha*x)
            x=(L_second*xopt/G_second)-W*Xw/G_second
            Xopt=Xopt+[xopt]
        else:
            N=len(Xopt)-1
    return print('Nmin=',Nmin,'     Rmin=',Rmin,'     N=',N,'     NFeed=',Nfeed,'     W=',W,'     D=',D)

Mc_Cabe(100, 0.35, 0.95, 0.05, 5, 2.5, 0)


""""""""""""""""" 
>--------------  Functions: 
    
           1|
            | enteshar ya nofoze dar gases >>.
            
             2| 
              | flux moli enteghal germ dar gases >>.
              
               3| 
                | shedat entegal germ >>. 
  
    
 '''The formulas that I used in these functions are taken from the Mass-Transfer-Operations-Robert-Treybal book.'''
--------------<  
  
>-------------- Convertors: 
 
    
           1-</
             /> Celsius_To_Kelvin (And vice versa)

            2-</
              /> CM**2/Centi_Second_to_M**2/Secend (And vice versa)
              
             3-</
               /> Atmosphere_to_Pascal (And vice versa)     
    
    
---------------<
    
""""""""""""""""" 






"""     1|
         | enteshar ya nofoze dar gases >> """

'''This function is used to calculate diffusivity in gases '''

# DAB = diffusivity, m2/s
# T = absolute temperature, K
# MAB,MAB = molecular weight of A and B, respectively, kg/kmol
# Pt = abs pressure, N/m2
# rAB = molecular separation at collision, nm = (r A + rB)/2
# εAB = energy of molecular attraction= (εAεB)**1/2
# k = Boltzmann's constant
#f(kT / εAB) collision function given by Fig. 2.5 in Mass-Transfer-Operations-Robert-Treybal

def Diffusivity(MA, MB, T, rAB, K, εAB, f, Pt=1.013*10**5):
    
    M = { 'Carbon' '(C)': 0.01201, #kg/kmol
 'Hydrogen' '(H)': 0.001008, #kg/kmol
 'Chlorine' '(Cl)': 0.03545, #kg/kmol
 'Bromine' '(Br)': 0.07990, #kg/kmol
 'Iodine' '(I)': 0.12690, #kg/kmol
 'Sulfur' '(S)': 0.03207, #kg/kmol
 'Nitrogen' '(N)': 0.01401, #kg/kmol
 'H2': 0.002016, #kg/kmol
 'O2': 0.03200, #kg/kmol
 'N2': 0.02802, #kg/kmol
 'Air': 28.97, #kg/kmol for dry air
 'CO': 0.02801, #kg/kmol
 'CO2': 0.04401, #kg/kmol
 'SO2': 0.06407, #kg/kmol
 'NO': 0.03001, #kg/kmol
 'N2O': 0.04402, #kg/kmol
 'Oxygen': 0.03200, #kg/kmol
 'NH3': 0.01703, #kg/kmol
 'H2O': 0.01802, #kg/kmol
 'H2S': 0.03408, #kg/kmol
 'COS': 0.06007, #kg/kmol
 'Cl2': 0.07090, #kg/kmol
 'Br2': 0.15980, #kg/kmol
 'I2': 0.25380, #kg/kmol
 'C3H6O':58.08  #kg/kmol
}

    if MA in M and MB in M:
        AB = float(((((((1/M[MA])**0.5 + (1/M[MB])**0.5) * 0.249) - 1.084) * (10**-4)) * T**(3/2)) * ((1/M[MA])**0.5 + (1/M[MB])**0.5)) 
        DAB = float(( AB / (((rAB**2) * Pt) * ((K*T / εAB) * f)))*-1)  #formula hame emtehan shode va ba data bedast amadeh be sorat mohasebeh dasti moghayeseh shodeh va ba deghat balayi javab doroost bodeh.
    
        return DAB
    else:
        return "Invalid inputs"
                 # Diffusivity (   MA  ,  MB  ,   T   ,  rAB ,   K  ,   εAB   ,  f   , Pt=1.013*10**5)
Diffusivity_value = Diffusivity('C3H6O', 'Air', 273.15, 0.495, 329.2, 58298.76, 0.389, 1.013*10**5)
print(Diffusivity_value,'M2/S') #inja mishe in khat va (Diffusivity_value) ro hazf kard chon tabe khorojii dareh va faghat mikhastam javabro ba vahedesh benevise:D







#=====================================================|
#=====================================================|
#=====================================================| 


"""      2| 
          | flux moli enteghal germ dar gases >> """ 

'''This function is used for calculation Molecular Diffusion in Gases, with two Type of Diffusion'''

#T=Kelvin
#Separate_Panels_Distance= Meter
#R=8.34 J/Mol.K 
#Pt=1.013 KG/M.S*2 , Pa , N/M*2 
#Diffusivity(DAB)= M**2/S 
#( Pressure_in_Panel_A, Pressure_in_Panel_B, PressureBM)= Pa or Atm 

def Diffusion_in_Gases(Diffusion_Type, Diffusivity, Separate_Panels_Distance, Pressure_in_Panel_A, Pressure_in_Panel_B, PressureBM, Gas_Constant, Total_Pressure=1.013*10**5):
   
                                                       
    DAB = {'H2-CH4': 6.25*10**-5, #M**2/S
           'O2-N2': 1.81*10**-5,                       
           'CO-O2': 1.85*10**-5,
           'CO2-O2': 1.39*10**-5,
           'Air-NH3': 1.98*10**-5,
           'Air-H2O': 2.58*10**-5,
           'Air-ethanol': 1.02*10**-5,
           'Air-ethyl-acetate': 0.87*10**-5,
           'Air-aniline': 0.74*10**-5,
           'Air-chlorobenzene': 0.74*10**-5,
           'Air-toluene': 0.86*10**-5}
    
    Temp = {'H2-CH4': 273.15,     #Kelvin
            'O2-N2': 273.15,
            'CO-O2': 273.15,
            'CO2-O2': 273.15,
            'Air-NH3': 273.15,
            'Air-H2O': 299.05,
            'Air-ethanol': 299.05,
            'Air-ethyl-acetate': 274.2,
            'Air-aniline': 299.05,
            'Air-chlorobenzene': 299.05,
            'Air-toluene': 299.05}
   
    # We can also write like this :
        
# if Diffusivity=='H2-CH4':
#        DAB=6.25*10**-5 
#       Temp=273.15
#if Diffusivity=='O2-N2':
#        DAB=1.81*10**-5 
#        Temp=273.15
#if Diffusivity=='CO-O2':
#        DAB=1.85*10**-5 
#        Temp=273.15
#if Diffusivity=='CO2-O2':
#        DAB=1.39*10**-5 
#        Temp=273.15
#if Diffusivity=='Air-NH3':
#        DAB=1.98*10**-5
#        Temp=273.15
#if Diffusivity=='Air-H2O':
#        DAB=2.58*10**-5 
#        Temp=299.05
#if Diffusivity=='Air-ethanol':
#        DAB=1.02*10**-5
#        Temp=299.05
#if Diffusivity=='Air-ethyl-acetate':
#        DAB=0.87*10**-5
#        Temp=274.2
#if Diffusivity=='Air-aniline':
#        DAB=0.74*10**-5
#        Temp=299.05
#if Diffusivity=='Air-chlorobenzene':
#        DAB=0.74*10**-5
#        Temp=299.05 
#if Diffusivity=='Air-toluene':
#        DAB=0.86*10**-5
#        Temp=299.05 

        
    if Diffusion_Type == 'Two_way_equimolal_diffusion':
        NA_Two = float(DAB[Diffusivity.replace('_', '-' )] / (Gas_Constant * Temp[Diffusivity.replace('_', '-' )] * Separate_Panels_Distance)) * (Pressure_in_Panel_A - Pressure_in_Panel_B) #( [Diffusivity.replace('_', '-' )] )-->>chon khodam eshtebahan az (_) bejaye (-) estephade karadam va error dad in bakhsh ro ezafe kardam ta agar bedalil tashaboh be eshtebah az (_) estefade shod, ba (-) jayegozin kone ta moshkel bartaraf she va edame bedeh. #Air_H2O--->>'Air-H2O'
        C = NA_Two
    elif Diffusion_Type == 'one-way_diffusion':
        NA_One = float((DAB[Diffusivity.replace( '_', '-' )] * Total_Pressure) / (Gas_Constant * Temp[Diffusivity.replace('_', '-' )] * Separate_Panels_Distance * PressureBM)) * (Pressure_in_Panel_A - Pressure_in_Panel_B)
        C = NA_One

    return C
                            #(  Diffusion_Type  , Diffusivity,(S_P_D),(P_A),(P_B),(PBM),Gas_Constant, Total_Pressure=1.013*10**5)
result_1 = Diffusion_in_Gases('one-way_diffusion', 'Air_H2O', 2*10**-3, 0.2, 0.1, 0.849, 8.314*10**3, 1.013*10**5)  #(Air_H2O--->>'Air-H2O') #((Pressure_in_Panel_A, Pressure_in_Panel_B, PressureBM) If these are written in Atmospheric units, there is no need to convert the unit to pa because it is simplified in the formula.). #(Because we want the (NA) unit to be based on ( KMOL/M**2 * S ), so we write the gas constant as (J/KMOL * K) which means we multiply it by 10**3)
print(result_1,'KMol/M**2 * S') #inja mishe in khat va (result_1) ro hazf kard chon tabe khorojii dareh va faghat mikhastam javabro ba vahedesh benevise:D .                                        




#=====================================================|
#=====================================================|
#=====================================================| 

"""       3| 
           |shedat entegal germ >> """ 
           
'''This function is used for Mass Transfer Intensity'''
           
#result_1=Diffusion_in_Gases(NA)
#The unit of (size_of_the_surface) is Meter

def Mass_Transfer_Intensity(result_1, size_of_the_surface): 
    MTI_ἠ = float(result_1 * size_of_the_surface)    #The size of the surface where mass transfer occurs(unite Meter).

    return MTI_ἠ

result_2 = Mass_Transfer_Intensity(result_1, 1) #inja andaze sathe entegal germ ro 1_square_meter dar nazar greftam
print(result_2,'KMol/s') #inja ham hamintor mishe in khat ro hazf kard a bejaye (result_1) bezarim (Diffusion_in_Gases).




#=====================================================|
#=====================================================|
#=====================================================|
"""
       / Convertors 
        
       / 1- Celsius_To_Kelvin (And vice versa)
        
       / 2- CM**2/Centi_Second__to__M**2/Secend (And vice versa)
        
       / 3- Atmosphere__to__Pascal (And vice versa)  
        
        """ 
        
''' 1-</
      /> Celsius_To_Kelvin (And vice versa) '''
      
'''This functions is use for converting Celsius to Kelvin and vice versa'''

def Celsius_To_Kelvin (C):
    K = float(C + 273.15) 
    return K 
Celsius_To_Kelvin (0)

#vice versa_______________

def Kelvin_To_Celsius (K): 
    C = float(K - 273.15)
    return C
Kelvin_To_Celsius (0)    

#=====================================================#    

'''2-</
     /> CM**2/Centi_Second_to_M**2/Secend (And vice versa)'''
     
'''This functions is use for convert square_centimeter/centi_second to square_meter/second and vice versa'''

def convert_cm2_cs_to_m2_s (cm2_cs) :
    m2_s = float(cm2_cs / 100)  # Conversion factor for area: 1 m^2 = 10000 cm^2  # Conversion factor for time: 1 s = 100 cs  # (10**-2) Combined conversion factor for square meters per second to square centimeters per centisecond
    return m2_s
convert_cm2_cs_to_m2_s(400)

#vice versa_________________________

def convert_m2_s_to_cm2_cs (m2_s) :
    cm2_cs = float(m2_s * 100)
    return cm2_cs
convert_m2_s_to_cm2_cs (4.0)
    
#=====================================================#     
                                                                                                           
'''3-</
     /> Atmosphere_to_Pascal (And vice versa) '''   
     
'''This functions is used for convert Atmosphere to Pascal and vice versa'''    

def Atmosphere_to_Pascal (atm):
    Pa= float(atm * 101325)
    return Pa 
Atmosphere_to_Pascal(1)
    
#vice versa___________________

def Pascal_to_Atmosphere (Pa):
    atm = float(Pa / 101325)
    return atm 
Pascal_to_Atmosphere (101325.0)
'''__This program is designed to calculate for composite materials
        includ : Black_Composite
                 White_Composite

   __Next step is to convert moles to percentages and percentages to moles.
   
'''

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

def Dental_Composite(totalR, totalF, totalC,
                     percentages_Filler, percentages_Resin, percentages_Colors,
                     molar_Si_SiO2, molar_O_SiO2,
                     molar_Ba_BaSiO3, molar_Si_BaSiO3, molar_O_BaSiO3,
                     molar_Zr_ZrO2, molar_O_ZrO2,
                     molar_Ti_TiO2, molar_O_TiO2
                     ):
    # Calculate the values corresponding to the percentage of each material
    percentages_dict = {}

    # Calculate the values of resin materials
    for i in range(len(Resin_Materials)):
        percentages_dict[Resin_Materials[i]] = totalR * (percentages_Resin[i] / 100)

    # Calculate the values of filler materials
    for i in range(len(Filler_Materials)):
        percentages_dict[Filler_Materials[i]] = totalF * (percentages_Filler[i] / 100)

    for i in range(len(Colors_Materials)):
        percentages_dict[Colors_Materials[i]] = totalC * (percentages_Colors[i] / 100)

    total_molar_SiO2 = molar_Si_SiO2 + 2 * molar_O_SiO2
    # Calculate the molar percentages of each element
    silicon_percentage = molar_Si_SiO2 / total_molar_SiO2
    oxygen_percentage = (2 * molar_O_SiO2) / total_molar_SiO2

    # Calculate the total molar mass of the filler material
    total_molar_BaSiO3 = molar_Ba_BaSiO3 + molar_Si_BaSiO3 + (3 * molar_O_BaSiO3)
    # Calculate the molar percentages of barium, silicon, and oxygen in the filler material
    barium_percentage = molar_Ba_BaSiO3 / total_molar_BaSiO3
    sil_Ba_percentage = molar_Si_BaSiO3 / total_molar_BaSiO3
    oxy_Ba_percentage = (3 * molar_O_BaSiO3) / total_molar_BaSiO3

    # Calculate the total molar mass of zirconium dioxide (ZrO2)
    total_molar_ZrO2 = molar_Zr_ZrO2 + (2 * molar_O_ZrO2)
    # Calculate the molar percentages of zirconium and oxygen in ZrO2
    zirconium_percentage = molar_Zr_ZrO2 / total_molar_ZrO2
    oxy_Zr_percentage = (2 * molar_O_ZrO2) / total_molar_ZrO2

    # Calculate the total molar mass of titanium dioxide (TiO2)
    total_molar_TiO2 = molar_Ti_TiO2 + (2 * molar_O_TiO2)
    # Calculate the molar percentages of titanium and oxygen in TiO2
    titanium_percentage = molar_Ti_TiO2 / total_molar_TiO2
    oxy_Ti_percentage = (2 * molar_O_TiO2) / total_molar_TiO2

    # Add calculated percentages to the dictionary
    percentages_dict['Silicon in SiO2'] = silicon_percentage
    percentages_dict['Oxygen in SiO2'] = oxygen_percentage
    percentages_dict['Barium in BaSiO3'] = barium_percentage
    percentages_dict['Silicon in BaSiO3'] = sil_Ba_percentage
    percentages_dict['Oxygen in BaSiO3'] = oxy_Ba_percentage
    percentages_dict['Zirconium in ZrO2'] = zirconium_percentage
    percentages_dict['Oxygen in ZrO2'] = oxy_Zr_percentage
    percentages_dict['Titanium in TiO2'] = titanium_percentage
    percentages_dict['Oxygen in TiO2'] = oxy_Ti_percentage

    return percentages_dict

# Values and percentages
totalR = 55
totalF = 35
totalC = 10
percentages_Filler = [30, 50, 20]
percentages_Resin = [65, 20, 15]
percentages_Colors = [70, 30]

# Calculate percentages
result = Dental_Composite(totalR, totalF, totalC,
                          percentages_Filler, percentages_Resin, percentages_Colors,
                          G_Mol_Si, G_Mol_O,
                          G_Mol_Ba, G_Mol_Si, G_Mol_O,
                          G_Mol_Zr, G_Mol_O,
                          G_Mol_Ti, G_Mol_O
                          )
print('\n\t\tStart to Dental_Composite Calcut\n ')
# Print the values
for material, value in result.items():
    print(f"Value_Percent --> {material}:", str(value) + "%" "\n")
    
print('\t\tComplited Dental_Composite Calcut ')    
print('\t\tStart convet to percentages_to_moles  \n')  

def percentages_to_moles(total, percentages):
    # Define the molecular weight of the composite mixture
    molar_weight = {'TEGDMA': 156.27, 'BIS_GMA': 512.67, 'UDMA': 398.48,
                    'Silica dioxide': 60.08, 'Barium silicate': 233.39, 'Zirconium dioxide': 123.22
                    }

    # Calculate the moles of each material
    moles = {}
    for material, percent in percentages.items():
        moles[material] = (percent / 100) * (total / molar_weight[material])

    return moles

# Values and percentages
total = 100
percentages = {'TEGDMA': 30, 'BIS_GMA': 50, 'UDMA': 20,
               'Silica dioxide': 65, 'Barium silicate': 20, 'Zirconium dioxide': 15
               }

# Convert percentages to moles
moles_result = percentages_to_moles(total, percentages)

# Print the moles of each material
for material, moles in moles_result.items():
    print(f"Moles of {material}:", moles)


print('\n\t\tComplited percentages_to_moles ')
print('\t\tStart convet to  moles_to_percentages  \n') 

def moles_to_percentages(total, moles):
    # Define the molecular weight of the composite mixture
    molar_weight = {'TEGDMA': 156.27, 'BIS_GMA': 512.67, 'UDMA': 398.48,
                    'Silica dioxide': 60.08, 'Barium silicate': 233.39, 'Zirconium dioxide': 123.22
                    }

    # Calculate the percentages of each material
    percentages = {}
    for material, mole in moles.items():
        percentages[material] = (mole * molar_weight[material] / total) * 100

    return percentages

# Values and moles
total_moles = {'TEGDMA': 0.1, 'BIS_GMA': 0.2, 'UDMA': 0.15,
               'Silica dioxide': 0.05, 'Barium silicate': 0.1, 'Zirconium dioxide': 0.2
               }

# Convert moles to percentages
percentages_result = moles_to_percentages(sum(total_moles.values()), total_moles)

# Print the percentages of each material
for material, percent in percentages_result.items():
    print(f"Percentages of {material}:", percent)
#variable
k= 534.6
pi= 3.14
F=96485
R=8.314
g=23


#function1
def Corrosion_Rate(num1,num2,num3,num4):
    '''
    

    Parameters
    ----------
    num1 : float
        the weight loss after exposure time(mg).
    num2 : float
        density(g/cm3).
    num3 : float
        exposed area(inch2).
    num4 : float
        time(hour).

    Returns 
    corrosin rate in mpy unit.
    

    '''
    Corrosion_Rate = (k*num1)/(num2*num3*num3)
    return Corrosion_Rate

Corrosion_Rate(2,3,4,)
#function2

def Contact_Angle(num1,num2,num3):
    '''
    

    Parameters
    ----------
    num1 : float
        a solid and atmospher surface tenssion.
    num2 : float
        a liquid and atmospher surface tension.
    num3 : float
        a solid and liquid boundary tension.

    Returns
    degrees

    '''
    import math
    costeta= (num1-num3)/num2
    teta = math.acos(costeta)
    teta_degree= teta*180/pi
    return teta_degree

Contact_Angle(12,10,4)
#function3
def Mass_Plating(num1,num2,num3,num4):
    '''
    

    Parameters
    ----------
    num1 : float
        current(A).
    num2 : float
        time(s).
    num3 : floar
        molar mass(g/mol).
    num4 : int
        number of electron.

    Returns
    precipitation mass in plating in g unit.

    '''
    m = num1*num2*num3/(num4*F)
    return m
Mass_Plating(1,4,3,2)
#convertor1
def Radians_to_degrees (num):
    '''
    
    This function is used for convert radians to degree
    '''
    degree=num*180/pi
    return degree
Radians_to_degrees(0.64)
#convertor2
def Weightpercent_to_ppm (num):
    '''
    
    This function is used for convert weight percent to ppm
    '''
    ppm= num*10000
    return ppm
Weightpercent_to_ppm(1)


    

#In the name of god
#A1_Haniyeh_salehi_project

'''
 baste be reshteye tahsilitoon: 
 parte 1:
 hade aghal 5 ta constant ( adade sabet ) tarif konid
'''
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

#--------------------------------

'''
part 2:
hade aghal 3 ta tabe tarif konid k b ESM, TOZIHAT, RETURN tavajo shavad
        
'''

def esm(vorodi1,vorodi2,vorodi32):
    '''
    
    **englishhhhhhh
    Parameters
    ----------
    vorodi1 : int
        a.
    vorodi2 : int
       mass.
    vorodi32 : TYPE
        DESCRIPTION.

    Returns --> khorojimon chie


    '''
    
    #body--------
    c=vorodi1+vorodi2
    
    
    return c #hatman return dashte bashe
#---------------------------
#FORMOL FESHAR
def Pressure(F,A):#dar in ja F,A ra be sorat voroodi minevisim
    P=F/A #taghsim mikonim
    return P# khoroji yek adad hast
 
J=Pressure(40,8) # in ja adad gozari anjam midahim
type(J)#Out[41]: float noee moshakhas mikonim
#FORMOL ENEGI POTANSIYEL GERANESHI
def G_P_E(m,g,h):#dar in ja vorodi migirim
    U=m*g*h #Amaliyat zarb anjam mishavad
    return U# khoroji yek adad hast
D=G_P_E(20,10,5) #in ja adad gozari anjam midahim

type(D)#Out[42]: int

#----------------------------------
'''
 part3:
 Dota convertor besaziid ke harkoodom baragshtesham besazid 
       ''' 
#2 ta tabe convertor ( harkodom rafto bargasht)

def M_to_mm (M):
    '''
    
    This function is usef for convert Meter to Milimeter
    '''
    mm=M*1000
    return mm

E=M_to_mm(15*1000)


#15000000

def mm_to_M (mm):
    '''
    
    This function is usef for convert Meter to Milimeter
    '''
    M=mm/1000
    return M
S=mm_to_M(15/1000)

#1.4999999999999999e-05
    
#------------------------------
#sakht do tbee convertor
#mohasebe tabdil fahrenheit be cantigerad
def Fahrenheit_to_C(F):
    '''
    This function is usef for convert Fahrenheit_to_Centigrade
    '''
    C=F-32/18 #fahrenheit menhaye32 mikonim baeed taghsim bar 18 ta cantigerad be dast biyad
    return C #dar in ja az return estefademikonim ta khoroji bedast ayad

T=Fahrenheit_to_C(45,32,18) # mohasebat ra anjam midahim va dar moteghayer R Zakhireh mikonim
#-----------------------------

def C_to_Fahrenheit(C):
    '''
    This function is usef for Centigrade_to_Fahrenheit
    '''
    F=C*1.8+32 #cantigerad ra zarbedar1.8 bealave32 mikonim ta  Fahrenheit be dast ayad
    return F #dar in ja az return estefademikonim ta khoroji bedast ayad
W=C_to_Fahrenheit(15,1.8,32) # mohasebat ra anjam midahim va dar moteghayer W Zakhireh mikonim

#---------------------------

        
Planck_constant=6.626*10**(-34) #js#
R=8.314 #J/mol.K#
Faraday_constant= 96485 #C/mol#
Boltzman_constany= 1.381*10**(-23) #J/K#
Avogardo_Num=6.022*10**(23) #mol^-1#



def Heat_Capacity (m,c,T1,T2):
    '''
    This function caclulate the amount of heat capacity of a mass
    Parameters
    ----------
    m : Float
        mass.
    c : float
        specific heat coefficient.
    T1 : float
        primary temperature.
    T2 : float
        secondary temperature.

    '''
    Q=m*c*(T2-T1)
    if Q<0:
        str='The mentioned material has lost heat'
    if Q==0:
        str='The mentioned material,s heat has not changed'
    if Q>0:
        str='The mentioned material has gained heat'
    return Q,str

##-----------------------------------------------------------------------------
import math
def Bragg_Law (h,k,l,a,y):
    '''
    This function calculate the diffraction angle of a incidence wavelenght through a crytal special plate

    Parameters
    ----------
    h : int
        x direction of the plate.
    k : int
        y direction of the plate.
    l : int
        z direction of the plate.
    a : float
        Unit Cell.
    y : float
       incidence wavelenght.

    '''
    square_sin_teta=(((y**2)*((h**2)+(k**2)+(l**2))))/4*a**2
    teta=math.asin((square_sin_teta)**1/2)
    
    return teta*180/math.pi
##-----------------------------------------------------------------------------

def Beer_Lambert_Law (a,l,I0):
    '''
    This function caclculate the output light intensity when a light pass through a material
    Parameters
    ----------
    a : float
        absorption coefficient.
    l : float
        length.
    I0 : float
        incomming light intensity.

    Returns
    -------
    I : float
        output light intensity.

    '''
    I=I0*(10**(-a*l))
    return I
    
##----------------------------------------------------------------------------

def Density_Convertor(d):
    '''
    This function change Kg/m3 to g/cm3

    Parameters
    ----------
    d : float
        density in Kg/m3.

    Returns
    -------
    d2 : float
        g/cm3.

    '''
    d2=d/1000
    return d2
##-----------------------------------------------------------------------------
  
def Pascal_to_mmHg(p):
    '''
    This function convert pascal to mmHg

    Parameters
    ----------
    p : float
        pressure (Pa).

    Returns
    -------
    None.

    '''
    mmHg=p/2
    return mmHg
##----------------------------------------------------------------------------


    




####constants
s=1###ehtemal tamam pishamadha
a=9.8###shetab dar harekat soghoot azad
Income_Tax_Perecent=0.09
Iran_Stock_Market_Volatility=0.05
Santiago_Bernabeu_Stadium_Crowd=78000


###
def Binomial_Probability(n,k,p):
    """
    

    Parameters
    ----------
    p :float
        ehtemal voghoo beyn 0ta 1
    n:int
        tedad dafat azmayesh
    K:int
        tedad dafaat bord

    return ehtemal rokh dadan pishamad dar fazaye do pishamadi(bord,shekast)
     ba tedad n azmayesh va tedad k  bord


    """
    q=1-p
    return p**k*q**(n-k)

#####
def Pythagorean_Relation(a,b,c):
    """
    

    Parameters
    ----------
    a : float
        tool zele aval.
    b : float
        tool zele dovom.
    c : float
        tool zele sevom.

    Returns be ma mige in adad mitoonan mosalas ghaem alazviye tashkil bedan 
    ba "mosalas ghaemalzaviye mishe" va" mosalas ghaemalzaviye nemishe"
    -------

    """
    if a**2==b**2+c**2 or b**2==a**2+c**2 or c**2==a**2+b**2:
        return "mosalas ghaemalzaviye mishe"
    else:
        return "mosalas ghaemalzaviye nemishe"

######
def Income_Tax(a):
    """
    

    Parameters
    ----------
    a : float
         daramad fard dar yek mah

    Returns maliyat bar daramad fard dar yek mah dar iran
    -------


    """
    return (a-10000000)*0.09


#####

#######
def Polygonal_Diameters(n):
    """
    

    Parameters
    ----------
    n : int
        tedad azlae chand zelei(az 2 ta bishtar).

    Returns tedad ghotrhaye chand zelei
    -------

    """
    return((n*(n-3))/2)
#######
def Velocity_Equation(V1,V2,a):
    """
    

    Parameters
    ----------
    V1 : float
        sorat avaliye moteharek.
    V2 : float
        sorat nahayi moteharek.
    a : floar
        shetab moteharek dar harekat ba shetab sabet.

    Returns mizan jabejayi dar harekat ba shetab sabet
    -------
   

    """
    return (V2**2-V1**2)/2*a
   
####convertor1
def Kmph_To_Mps(V1):
    """
    This function is uesd to convert Kilometer per  hour to meter per second
     """
    
    V2=V1/3.6
    return V2

####reverse
def Mps_To_Kmph(V1):
    """
    This function is used to convert meter per second to kilometer per hour
    """
    V2=V1*3.6
    return V2
#####convertor2
def Pascal_To_CmHg(P1):
    """
    This function is used to convert Pascal   to centimeter mercury 
    """
    P2=P1/1333.22
    return P2
#####reverse
def CmHg_To_Pascal(P1):
    """
    This function is used to convert mercury centimeter to Pascal
    """
    P2=P1*1333.22
    return P2
    
    
    
"""
A1_Shahrokh_Akbari

send to ai.2024.pilehvar@gmail.com

@author: Shahrokh Akbari sis

Modelling of Carbon Diffusion Process _ Heat Treatment
Function 1 : diffusion Co-efficient 
Function 2 : Van-Ostrand-Dewey solution to fick's second Law
(using boundary conditions of infinite source and length but finite maximum carbon diffusivity)
Semi-infinite condition
Function 3 : constant amount of carbon in the system (thin source)
"""
#The constants:
#sample size:
a=7
b=10
c=60
#Diffusion Coeffitent
Do=0.12
#temperature in kelvin
T=1293
#different measures for time used as example in second[s]
ti=3600
tm=10800
tf=18000
#carbon diffusivity at Tf cm**2/s
Df=1.97*10**(0-8)
#surface carbon content is Cs
Cs=0.16
Co=0.4
#for formula :
"""
-E/r = 16000
"""
#average activation Energy of carbon in steel
Er=16000
Ec=20200
#Gas constant
R=8.314
#------------------------------------------------------------------------------
#Function 1
def Diffusion_Coefficient(T):
    """
    Parameters
    ----------
    T : integer
        the Temperature for which we require the Coeffiecient.
    Returns
    -------
    Dm : integer
        Diffusion Coefffecient at the required Temperature.

    """
    global Dm
    Dm=Do*2.718**(-Ec/T)
    return Dm
#the data tested in the Journal and the Research paper equates the data generated by the Function 2
print(Diffusion_Coefficient(T))
#------------------------------------------------------------------------------
#Function 2
"""
here i had to use the import math formula because i needed the erf() function.
there exists a function for erf() in the mentioned excel sheet
"""
#for this function we need Dm from the first function i did not add it seperately but it is possible to add a new input for the function taking the diffuion coefficient from us.
def Carbon_Content(x,t):
    """
    Parameters
    ----------
    x : integer
        Distance from surface (cm)
    t : integer
        Time given for diffusion
    Returns
    -------
    Cxt : integer
        Carbon content at the required time and distance from surface [thickness]

    """
    import math
    Cxt=Cs+(Co-Cs)*math.erf(x/2/Dm**0.5/t**0.5)
    return Cxt
print(Carbon_Content(10, ti))
#------------------------------------------------------------------------------
#Function 3
def Carbon_content(B,x,t):
    """
    Parameters
    ----------
    B : integer
        the total available carbon content
    x : integer
        the thickness in question
    t : ineteger
        time 
    Returns
    -------
    Cxt : integer
        carbon content at the required time and location according to total carbon provided

    """
    Cxt=(B/(3.14*Dm*t)**0.5)*2.718**(-(x**2)/4/Dm/t)
    return Cxt
print(Carbon_content(0.019,0,ti))
#B is a constant that is usually between 1 and 2 depoendant on the type of steel used
#------------------------------------------------------------------------------
#convertor1
def sec_to_hour(t):
    t=t/3600
    return t
print(sec_to_hour(10800))
#convertor 1 reversed
def hour_to_sec(t):
    t=t*3600
    return t
print(hour_to_sec(5))
#------------------------------------------------------------------------------
#convertor 2 
def Cm_to_um(x):
    a=x*10000
    return a
print(Cm_to_um(0.002))
#convertor 2 reverse
def um_to_Cm(x):
    a=x/10000
    return a
print(um_to_Cm(500))




#===============================
#1_Constants:

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





#===========================================
#2_Functions:

    
#First_Function------>Incompressible_Fluids_Pressure
'''
This function gets the density of the fluid and the depth in fluid we are calculating the pressure for as inputs,
and returns the local pressure of that incmpressible fluid at that specofoc depth as the output.
The units for variables are:
P and P_0------>Pa
d------>Kg/m^3
h------>m
'''
def Incompressible_Fluids_Pressure(d,h):  
    P_0=101325#(Pa)
    g=9.81#(m/s^2)
    P=d*g*h+P_0
    return P





#Second_Functin------>Bouyancy_Force
'''
This function calculates the bouyancy force that is applied to an object by a fluid.
The input d is the density of the fluid with unit of Kg/m^3 
and V is the volume of that part of the object which has moved or you can say is inside of the fluid with the unit of m^3
The functon returns F_b(bouyancy force) with the unit of N or in SI form--->Kg.m/s^2
'''
def Bouyancy_Force(d,V):
    g=9.81#(m/s^2)
    F_b=d*g*V
    return F_b






#Third_Function------>Reynolds_Number
'''
Reynods number is a dimensionless value that can help us to find out if a fluid flow is torbulent or laminar.
This function gives the Re(Reynolds Number) for a fluid flow in a pipe as an output. 
In order to calculate the value of Re we need these variables as inputs:
d as density of fluid[Kg/m^3]
u as flow speed[m/s]
vis as viscosity of the fluid[Pa.s or SI:Kg/m.s]
D_H as Hydraulic Diameter*(characterisitic length) [m]

*hydraulic diameter is the inside diameter of the pipe.

!!!note that this function is only for cylenders and in any other container with a cross section other than a circle
you should define a new characteristic length.
'''   
def Reynolds_Number_Pipe(d,u,vis,D_H):
    Re=d*u*D_H/vis
    return Re




#Convertors:

#First_convertor------>Electronvolt_to_Joule
'''
This function converts the values of energy from electronvolt to joule.
'''
def Electronvolt_to_Joule(e_v):
    Joule=e_v*1.6022e-19
    return Joule


#Second_convertor------>Degree_to_Radian
'''
This function converts values of angle from degree to radian.
'''
def Degree_to_Radian(deg):
    rad=deg*3.141592653589793/180
    return rad









#constants
Latent_heat = 1.16 * (10 ** 9)
Gama_sl = 0.132
T_m = 1064 #in centigrade
A = 0.413 #in nano meter
Pi = 3.14
K = 1.38 * (10 ** (-23)) #Boltzmann constant
N_t = 3 * (10 ** 33)


#convertors
def Meter_to_nanometer(value):
    return value * (10**9)
    
def Nanometer_to_meter(value):
    return value * (10 ** (-9))

def Centigrade_to_kelvin(value):
    return value + 237

#functions
def Critical_radius(Temp_difference):
   return (2 * Gama_sl * T_m)/(Latent_heat * Temp_difference)

def Activation_energy(Temp_difference):
    return (16 * Pi * (Gama_sl ** 3) * (T_m ** 2))/(3 * (Latent_heat ** 2) * (Temp_difference ** 2))

def Number_of_unit_cells(Radius):
    return (4 * Pi * (Radius ** 3))/(3 * (converted_A ** 3))

def Nucleation(Temprature, activation):
    return (N_t * math.exp(-(activation/(K * Temprature))))

#inputs
Initial_temp = int(input("Enter temprature in centigrade: "))
Temp_difference = T_m - Initial_temp

#calculations
converted_A = Nanometer_to_meter(A)
converted_Initial_temp = Centigrade_to_kelvin(Initial_temp)

Gold_critical_radius = Critical_radius(Temp_difference)
Gold_activation_energy = Activation_energy(Temp_difference)
Gold_unit_cells = Number_of_unit_cells(Gold_critical_radius)
Gold_nucleation = Nucleation(converted_Initial_temp, Gold_activation_energy)


#Printed out values
print("This is gold's critical radius: ",Gold_critical_radius)
print("This is gold's critical activation energy: ",Gold_activation_energy)
print("This is gold's unit cells: ",Gold_unit_cells)
print("Number of gold nucleations per cubic meter : ",Gold_nucleation)


N_a = 6.022*10**23              #molecules/mol (Avogadro's Number)
k = 1.381*10**(-23)             #J/K (Boltzmann's Constant)
alpha = 5.67*10**(-8)           #W/m^2/K^4 (Stefan-Boltzmann Constant)
R = 8.314                       #J/mol/K (Gas costant)
Visible = range(380, 750, 10)   #Visible Wavelength
h = 6.626*10**(-34)             #J.s (Planck's Constant)
c = 2.998*10**8                 #m/s (speed of light)

def Copolymer_Type(Copolymer, Polymer_num=2):
    '''
    Copolymer: str
    Unit types are shown by charachters which are seperated by a space
    Example: 'polymer_A polymer_A polymer_B polymer_A polymer_A polymer_B'
    
    Polymer_num: int
    It represents the number of polymers in this structure which could be 2 or 3
    default = 2
    '''
    valid = [2, 3]
    if Polymer_num not in valid:
        raise ValueError('Polymer_num should be 2 or 3')
    
    co_list = Copolymer.split()
    change = []
    unique = 0
    unique_list = []
    for i in range(1,len(co_list)):
        if co_list[i-1] not in unique_list:
            unique_list.append(co_list[i-1])
            unique += 1
        if co_list[i] != co_list[i-1]:
            change.append(True)
        else:
            change.append(False)
            
    if change.count(True) == 0:
        print('Error: It is not a copolymer') 
    elif Polymer_num != unique:
        print('Error: The number of unique units in the copolymer is not compatible with the Polymer_num entered')
    elif Polymer_num == 2:
        if change.count(False) == 0:
            copolymer_type = 'Alternative'
            print(copolymer_type)
        elif change.count(True)==1:
            copolymer_type = 'Block'
            print(copolymer_type)
        else:
            copolymer_type = 'Random'
            print(copolymer_type)
        return copolymer_type
    
    


import numpy as np
    
def PengRobinson(T = None,P = None,Tc = None,Pc = None,w = None,MW = None,Phases = None):
    """
    PengRobinson.m : calculates the compressibility factor,fugacity coefficient and density
    of a pure compound with the Peng Robinson equation of state (PR EOS)

    Parameters
    ----------
    T : float
        Temperature [=] K
    P : float
        Presure [=] Pa
    Tc : float
        Critical temperature [=] K
    Pc : float
        Critical presure [=] Pa
    w : float
        Accentic factor
    MW : float
        Molar weigth [=] kg/mol.
    Phases : int
        if Phases == 1, then calculates liquid fugacity;
        if Phases == 0 then calculates vapor fugacity

    Returns
    -------
    Z : flout
        Compressibility factor
    fhi : float
        Fugacity coefficient
    density : float
        Density

    """
    R = 8.314
    
    # Reduced variables
    Tr = T / Tc
    
    # Parameters of the EOS for a pure component
    m = 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    alfa = (1 + m * (1 - np.sqrt(Tr))) ** 2
    a = 0.45724 * (R * Tc) ** 2 / Pc * alfa
    b = 0.0778 * R * Tc / Pc
    A = a * P / (R * T) ** 2
    B = b * P / (R * T)
    # Compressibility factor
    Z = np.roots(np.array([1,- (1 - B),(A - 3 * B ** 2 - 2 * B),- (A * B - B ** 2 - B ** 3)]))
    ZR = []
    for i in range(3):
        if type(Z[i])!='complex':
            ZR.append(Z[i])
    
    if Phases == 1:
        Z = np.amin(ZR)
    else:
        Z = np.amax(ZR)
    
    # Fugacity coefficient
    fhi = np.exp(Z - 1 - np.log(Z - B) - A / (2 * B * np.sqrt(2)) * np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)))
    if True:
        density = P * MW / (Z * R * T)
        result = np.array([Z,fhi,density])
    else:
        'No real solution for "fhi" is available in this phase'
        result = np.array(['N/A','N/A','N/A'])
    

    return Z,fhi,density



#part1
water_density=1000 #kg/m3
gas_constant=8.314*10**-5 #m3.bar/k.mol
Nacl_MW=58 #g/mol
PI=3.141592653589793238462643383279502884197
Material=['0=Aceton','1=Acetonitrile','2=Acrylonitrile','3=Ammonia','4=Aniline','5=Benzalehyde','6=Benzene','7=n-Butane','8=n-Butanol','9=iso-Butane','10=iso-Butanol','11=Butylacetate','12=Carbondisulphide','13=Carbontetrachloride','14=Chlorobenzene','15=Chloroform','16=Cyclohexane','17=Cyclohexanol','18=Cyclohexanone','19=Cyclopentane','20=Dioxane','21=Dichloromethane','22=Diethylether','23=Diethylamine','24=Ethanol','25=Ethylacetate','26=Ethylbenzene','27=Ethylamine','28=Formicacid','29=Furfural','30=n-Hexane','31=n-Heptane','32=Methanol','33=Methylacetate','34=Nitrobenzene','35=Nitrogen','36=n-Octane','37=Oxygen','38=Octanol','39=n-Pentane','40=Phenol','41=n-Propanol','42=iso_Propanol','43=Propane','44=Pyridine','45=Styrene','46=Tetrahydrofuran','47=Toluene','48=Trichloroethylene','49=Triethylamine','50=o-Xylene','51=p-Xylene','52=Water']
A=[16.39112,16.90395,15.92847,17.51202,16.67784,6.73163,15.9037,15.68151,17.62995,15.77506,18.02933,16.4145,15.77889,15.8434,16.4,16.017,15.7794,19.23534,16.40517,15.8602,17.1151,17.0635,16.5414,15.73382,18.68233,16.35578,16.04305,7.3862,15.9938,15.14517,15.9155,15.877,18.61042,16.58646,16.42172,15.3673,15.9635,15.06244,7.18653,15.8365,15.9614,17.8349,20.4463,15.7277,16.152,15.94618,16.11023,16.00531,15.01158,15.7212,7.00154,6.99052,18.5882]
B=[2787.5,3413.1,2782.21,2363.24,3858.22,1369.46,2789.01,2154.9,3367.12,2133.24,3413.34,3293.66,2585.12,2790.78,3485.35,2696.25,2778,5200.53,3677.63,2589.2,3579.78,3053.08,2847.72,2434.73,3667.7,2866.6,3291.66,1137.3,2982.45,2760.09,2738.42,2911.32,3392.57,2839.21,3485.35,648.59,3128.75,674.59,1515.427,2477.07,3183.67,3310.4,4628.95,1872.82,3124.45,3270.26,2768.37,3090.78,2345.48,2674.7,1476.393,1453.43,3984.92]
C=[229.67,250.48,222,250.54,200,177.081,220.79,238.74,188.7,245,199.97,210.75,236.46,226.46,224.87,226.24,223.14,251.7,212.7,231.36,240.35,252.6,253,212,226.1,217.9,213.8,235.85,218,162.8,226.2,226.65,230,228,224.84,270.02,209.85,263.07,156.767,233.21,159.5,198.5,252.64,250,212.66,206,226.3,219.14,192.73,205,213.872,215.307,233.43]

           
#part3

def Flux_convertor1(a):
    '''
This function is used for convert ml/cm2.min to L/m2hr
    a : flow in ml/cm2.s unit
    Returns flow in L/m2hr unit
    '''
    b=a*600
    return b
Flux_convertor1(113.3)
def Flux_convertor2(a):
    '''
This function is used for convert L/m2hr to ml/cm2.min
    a : flow in L/m2hr unit
    Returns flow in ml/cm2.min unit
    '''
    b=a/600
    return b
Flux_convertor2(1)

def Pressure_convertor1(a):
    '''
This function is used for convert mmHg to bar
    a : Pressure in mmHg unit
    Returns Pressure in bar unit
    '''
    b=a/760
    return b
Pressure_convertor1(760)
def Pressure_convertor2(a):
    '''
This function is used for convert bar to mmHg
    a : Pressure in bar unit
    Returns Pressure in mmHg unit
    '''
    b=a*760
    return b
Pressure_convertor2(1)

def Temperature_convertor1(a):
    '''
This function is used for convert Celsius to Fahrenheit
    a : Temperature in Celsius unit
    Returns Temperature in Fahrenheit unit
    '''
    b=(a*1.8)+32
    return b
Temperature_convertor1(10)
def Temperature_convertor2(a):
    '''
This function is used for convert Fahrenheit to Celsius
    a : Temperature in Fahrenheit unit
    Returns Temperature in Celsius unit
    '''
    b=(a-32)/(1.8)
    return b
Temperature_convertor2(10)

def Velocity_convertor1(a):
    '''
This function is used for convert m/s to km/hr
    a : Velocity in m/s unit
    Returns Velocity in km/hr unit
    '''
    b=a*3.6
    return b
Velocity_convertor1(10)
def Velocity_convertor2(a):
    '''
This function is used for convert km/hr to m/s
    a : Velocity in km/hr unit
    Returns Velocity in m/s unit
    '''
    b=a/3.6
    return b
Velocity_convertor2(10)


#part2

def Perm(r,F,P):
    '''
    This function is used for membrane permeability calculation 
    r : membrane radius (cm)
    F : Flowrate (ml/min)
    P : applied pressure (Bar)
    Returns membrane permeability in LMH/Bar Unit 
    '''
    Area=PI*(r**2)
    Flux=F/Area
    flux=Flux_convertor1(Flux)
    Perm=flux/P
    return Perm
permeability=Perm(3, 424.6, 20)

def Rejection(CF,CP):
    '''
    This function is used for membrane Rejection calculation 
    CF : Pollutant concentration in Feed (g/l)
    CP : Pollutant concentration in Permeate (g/l)
    Returns membrane Rejection (%)  
    '''
    Rej=1-(CP/CF)
    return Rej
rejection=Rejection(2.17, 1.96)

def OP(M,T):
    '''
    This function is used for membrane Osmotic Pressure calculation 
    M : NaCl Concentration in Feed (g/l)
    T : Temperature (c)
    Returns membrane Osmotic Pressure in Bar
    '''
    m=M/(Nacl_MW)
    t=T+273.15
    Osmotic_Pressure=1.8*gas_constant*0.01*m*t
    return Osmotic_Pressure
op=OP(5, 25)

def TB():
    '''
    This function is used for bubble temp calculation in mixed solution
    P : pressure (mmhg)
    N : number of component
    Returns bubble temp
    '''
    Material=['0=Aceton','1=Acetonitrile','2=Acrylonitrile','3=Ammonia','4=Aniline','5=Benzalehyde','6=Benzene','7=n-Butane','8=n-Butanol','9=iso-Butane','10=iso-Butanol','11=Butylacetate','12=Carbondisulphide','13=Carbontetrachloride','14=Chlorobenzene','15=Chloroform','16=Cyclohexane','17=Cyclohexanol','18=Cyclohexanone','19=Cyclopentane','20=Dioxane','21=Dichloromethane','22=Diethylether','23=Diethylamine','24=Ethanol','25=Ethylacetate','26=Ethylbenzene','27=Ethylamine','28=Formicacid','29=Furfural','30=n-Hexane','31=n-Heptane','32=Methanol','33=Methylacetate','34=Nitrobenzene','35=Nitrogen','36=n-Octane','37=Oxygen','38=Octanol','39=n-Pentane','40=Phenol','41=n-Propanol','42=iso_Propanol','43=Propane','44=Pyridine','45=Styrene','46=Tetrahydrofuran','47=Toluene','48=Trichloroethylene','49=Triethylamine','50=o-Xylene','51=p-Xylene','52=Water']
    A=[16.39112,16.90395,15.92847,17.51202,16.67784,6.73163,15.9037,15.68151,17.62995,15.77506,18.02933,16.4145,15.77889,15.8434,16.4,16.017,15.7794,19.23534,16.40517,15.8602,17.1151,17.0635,16.5414,15.73382,18.68233,16.35578,16.04305,7.3862,15.9938,15.14517,15.9155,15.877,18.61042,16.58646,16.42172,15.3673,15.9635,15.06244,7.18653,15.8365,15.9614,17.8349,20.4463,15.7277,16.152,15.94618,16.11023,16.00531,15.01158,15.7212,7.00154,6.99052,18.5882]
    B=[2787.5,3413.1,2782.21,2363.24,3858.22,1369.46,2789.01,2154.9,3367.12,2133.24,3413.34,3293.66,2585.12,2790.78,3485.35,2696.25,2778,5200.53,3677.63,2589.2,3579.78,3053.08,2847.72,2434.73,3667.7,2866.6,3291.66,1137.3,2982.45,2760.09,2738.42,2911.32,3392.57,2839.21,3485.35,648.59,3128.75,674.59,1515.427,2477.07,3183.67,3310.4,4628.95,1872.82,3124.45,3270.26,2768.37,3090.78,2345.48,2674.7,1476.393,1453.43,3984.92]
    C=[229.67,250.48,222,250.54,200,177.081,220.79,238.74,188.7,245,199.97,210.75,236.46,226.46,224.87,226.24,223.14,251.7,212.7,231.36,240.35,252.6,253,212,226.1,217.9,213.8,235.85,218,162.8,226.2,226.65,230,228,224.84,270.02,209.85,263.07,156.767,233.21,159.5,198.5,252.64,250,212.66,206,226.3,219.14,192.73,205,213.872,215.307,233.43]
    P=float(input('Pressure='))
    N=int(input('Number of component='))
    x=[]
    T=[]
    a=[]
    b=[]
    c=[]
    d=len(Material)
    Index=[]
    guess=0
    p=[]
    U=0
    import math
    print(Material)
    for i in range(0,N):
        I=[int(input('enter component index='))] 
        Index=Index+I
    i=0
    for e in range(0,d):
        if i==N:
            break
        elif Index[i]==e:
            a=a+[[A[e]]]
            b=b+[[B[e]]]
            c=c+[[C[e]]]
            i=i+1
        else:
            e=e+1
    for i in range(0,N):
        X=[float(input('enter component mole fraction='))] 
        x=x+X
    for i in range(0,N):
        t=[(B[Index[i]]/((A[Index[i]])-math.log(P,math.e))-C[Index[i]])]
        T=T+t
    n=0
    while n<100000:
        if int(U)==int(P):
            TBP=guess
        else:
            guess=guess+0.01
            p=[]
            for i in range(0,N):        
                pr=[math.exp(A[Index[i]]-((B[Index[i]])/((C[Index[i]])+guess)))]
                p=p+pr
                U=0
            for i in range(0,N):           
                u=p[i]*x[i]
                U=u+U
        n=n+1
    return TBP
BubbleTemp=TB()     

def TD():
    '''
    This function is used for dew temp calculation in mixed solution
    P : pressure (mmhg)
    N : number of component
    Returns bubble temp
    '''
    Material=['0=Aceton','1=Acetonitrile','2=Acrylonitrile','3=Ammonia','4=Aniline','5=Benzalehyde','6=Benzene','7=n-Butane','8=n-Butanol','9=iso-Butane','10=iso-Butanol','11=Butylacetate','12=Carbondisulphide','13=Carbontetrachloride','14=Chlorobenzene','15=Chloroform','16=Cyclohexane','17=Cyclohexanol','18=Cyclohexanone','19=Cyclopentane','20=Dioxane','21=Dichloromethane','22=Diethylether','23=Diethylamine','24=Ethanol','25=Ethylacetate','26=Ethylbenzene','27=Ethylamine','28=Formicacid','29=Furfural','30=n-Hexane','31=n-Heptane','32=Methanol','33=Methylacetate','34=Nitrobenzene','35=Nitrogen','36=n-Octane','37=Oxygen','38=Octanol','39=n-Pentane','40=Phenol','41=n-Propanol','42=iso_Propanol','43=Propane','44=Pyridine','45=Styrene','46=Tetrahydrofuran','47=Toluene','48=Trichloroethylene','49=Triethylamine','50=o-Xylene','51=p-Xylene','52=Water']
    A=[16.39112,16.90395,15.92847,17.51202,16.67784,6.73163,15.9037,15.68151,17.62995,15.77506,18.02933,16.4145,15.77889,15.8434,16.4,16.017,15.7794,19.23534,16.40517,15.8602,17.1151,17.0635,16.5414,15.73382,18.68233,16.35578,16.04305,7.3862,15.9938,15.14517,15.9155,15.877,18.61042,16.58646,16.42172,15.3673,15.9635,15.06244,7.18653,15.8365,15.9614,17.8349,20.4463,15.7277,16.152,15.94618,16.11023,16.00531,15.01158,15.7212,7.00154,6.99052,18.5882]
    B=[2787.5,3413.1,2782.21,2363.24,3858.22,1369.46,2789.01,2154.9,3367.12,2133.24,3413.34,3293.66,2585.12,2790.78,3485.35,2696.25,2778,5200.53,3677.63,2589.2,3579.78,3053.08,2847.72,2434.73,3667.7,2866.6,3291.66,1137.3,2982.45,2760.09,2738.42,2911.32,3392.57,2839.21,3485.35,648.59,3128.75,674.59,1515.427,2477.07,3183.67,3310.4,4628.95,1872.82,3124.45,3270.26,2768.37,3090.78,2345.48,2674.7,1476.393,1453.43,3984.92]
    C=[229.67,250.48,222,250.54,200,177.081,220.79,238.74,188.7,245,199.97,210.75,236.46,226.46,224.87,226.24,223.14,251.7,212.7,231.36,240.35,252.6,253,212,226.1,217.9,213.8,235.85,218,162.8,226.2,226.65,230,228,224.84,270.02,209.85,263.07,156.767,233.21,159.5,198.5,252.64,250,212.66,206,226.3,219.14,192.73,205,213.872,215.307,233.43]
    P=float(input('Pressure='))
    N=int(input('Number of component='))
    x=[]
    T=[]
    a=[]
    b=[]
    c=[]
    d=len(Material)
    Index=[]
    guess=0
    p=[]
    U=0.1
    import math
    print(Material)
    for i in range(0,N):
        I=[int(input('enter component index='))] 
        Index=Index+I
    i=0
    for e in range(0,d):
        if i==N:
            break
        elif Index[i]==e:
            a=a+[[A[e]]]
            b=b+[[B[e]]]
            c=c+[[C[e]]]
            i=i+1
        else:
            e=e+1
    for i in range(0,N):
        X=[float(input('enter component mole fraction='))] 
        x=x+X
    for i in range(0,N):
        t=[(B[Index[i]]/((A[Index[i]])-math.log(P,math.e))-C[Index[i]])]
        T=T+t
    n=0
    while n<100000:
        if int(1/U)==int(P):
            TDP=guess
        else:
            guess=guess+0.01
            p=[]
            for i in range(0,N):        
                pr=[math.exp(A[Index[i]]-((B[Index[i]])/((C[Index[i]])+guess)))]
                p=p+pr
                U=0
            for i in range(0,N):           
                u=x[i]/p[i]
                U=u+U
        n=n+1
    return TDP
DewTemp=TD() 

"""
Created on Thu Mar  7 23:48:24 2024

@author: Mojtaba Bakhshi Maldar

Discription:
    This field includes a wide range of skills that are mainly used in the field of planning and control. In fact, industrial engineering is a junction between engineering and management.
    One of the most used and practical areas of this field is project management and control. Therefore, this project is dedicated to creating a code to directly calculate the indicators of the Earned Value Method (EVM is a quantitative technique to calculate the indicators of project progress).
    this project contains 3 main parts:
        1.1."Industrial Engineering constanst" that maily used for calculating indicators.
        1.2."EVM functions" that give constants and return the values of indicators
        1.3."Convertors functions" for convert some units 

"""

#Part1 ------------------------------------------------------------------------
#IE_Constants

k=1.96              #k=Z(0.05) Normal Standard Distribution
pf=0.6209           #f(P/F , i% , n) = 0.6209 for i=10% and n=5
ag=3.1936           #f(A/G , i% , n) = 3.1936 for i=18% and n=10
t_stu=0.9277        #T-Student Distribution for x=1.5 and n=30 (degree freedom 29)
d=1.128             #Constant coefficients for n=2


#Part2 ------------------------------------------------------------------------
#EVM_Functions

'''
There are 4 main inputs in EVM:
    BAC Stands for Budget at Complete
    BCWS Known as PV stands for Budgeted Cost of Work Scheduled
    BCWP Known as EV stands for Budgeted Cost of Work Performed
    ACWP Known as AC stands for Actual Cost of Work Performed
Also, there are 2 main indicators which each of them has 2 outputs:
    Schdule Indicators: represent for Time criteria of project include:
        1.SV=EV-PV
        2.SPI=EV/PV
    Cost Indicators: represent for Cost criteria of project include:
        1.CV=EV-AC
        2.CPI=EV/AC
It should be mentioned that there are many other indicators (such as EAC) which are calculated according to the need.
'''

def Schedule_Indicators(pv,ev):
    
    global sv,spi
    sv=ev-pv
    spi=ev/pv
    
    return sv,spi


Schedule_Indicators(1000, 1200)   #Out[5]: (200, 1.2) SV=1200-1000=200 and SPI=1200/1000=1.2

def Cost_Indicators(ac,ev):
    
    global cv,cpi
    cv=ev-ac
    cpi=ev/ac
    
    return cv,cpi

Cost_Indicators(1500, 1200)     #Out[7]: (-300, 0.8) CV=1200-1500 and CPI=1200/1500=0.8

def Estimate(bac,ev,ac,pv):
    
    Cost_Indicators(ac, ev)
    Schedule_Indicators(pv, ev)
    pf=cpi*spi
    etc=(bac-ev)/pf
    eac=round(etc+ac)
    
    return eac

Estimate(10000, 1200, 1500, 1000)   #Out[23]: 10667 This value shows us how a project budget can change due to its performance during the life of the project.


#Part3 ------------------------------------------------------------------------
#Convertors

'''
In this part, we define functions that play the role of unit converters:
    1.Distance convertors like foot (ft) to mile (mi) and its reverse
    2.Data convertors like byte (b) to kilobyte (kb) and its reverse
'''

def Foot_to_Mile (ft):
    
    mi=0.000189393939*ft
    
    return mi

def Mile_to_Foot (mi):
    
    ft=5280*mi
    
    return ft

Foot_to_Mile(1000)  #Out[26]: 0.189393939
Mile_to_Foot(0.5)   #Out[27]: 2640.0

def Byte_to_Kilobyte (b):
    
    kb=0.0009765625*b
    
    return kb

def Kilobyte_to_Byte (kb):
    
    b=1024*kb
    
    return b

Byte_to_Kilobyte(1024)  #Out[29]: 1.0
Kilobyte_to_Byte(0.5)   #Out[30]: 512.0


###----------Part1---------------
R=8.314 ##gas constant in J/mol.K
D=2.3*(10**(-5)) ##oxygen diffusion coefficient in water μm2/s
S=28.34 ##standard enthropy for solid alumminuim in J/mol.K
k=1.38*(10**(-23))  ##boltzmann constant for gas particle kinetic energy
Tm=14.025 ## hydrojen melting point in K
###---------------Part2-------------------------------------------------------
##part2-section1*********
##polymer burning rate calculation function
def Burning_Rate(L,t):##L=burning length in mm & t=burning time in sec
    V=(60*L)/t ##V= burning rate in mm/sec according to the ASTM D3801 UL-94 test
    return V

Burning_Rate(500,190)

##part2--section2***********
##polymer crystallinity calculation in DSC analysis
def Crystal_Percent(H,W,H100):#H=polymer enthalpy in mJ & W=polymer weight in mg & H100=neede enthalpy for polymer with 100%crystallinity in j/gr
    Xc=((H/W)/H100)*100 #Xc=polymer crystallinity in %
    return Xc

Crystal_Percent(816.2,11.8,203)##calculation for Neat eva

##Part2-section3***********
##flame retardant wt% in polymer matrix
def Filler_Weight(M,FR1,FR2,FR3):#M=polymer matrix weight in gr & F1=first flame retardant weight in gr & F2=second flame retardant weight in gr & F3=third flame retardant weight in gr
    a=[FR1/(FR1+FR2+FR3+M)]#FR1 weight%
    b=[FR2/(FR1+FR2+FR3+M)]#FR2 weight %
    c=[FR3/(FR1+FR2+FR3+M)]#FR3 weight %
    return a,b,c

Filler_Weight(25,1,10,5)

##-----------------------------------------part3-------------------------------
##part3--section1*************************
##converting celsius temperature to farenheit
def Celesius_to_Farenheit(a):##a= temperature in celsius
    b=(1.8*a)+32##b=temperature in farenheit
    return b

Celesius_to_Farenheit(50)
##convertingfarenheit temperature to celsius
def Farenheit_to_Celsius(a):##a=temperature in farenheit
    b=(a-32)/1.8##b=temperature in celsius
    return b

Farenheit_to_Celsius(122)

##part3---section2*********************
##converting ppm to percent
def Ppm_to_Percent(a):##a=ion concentration in ppm in brine
    b=a/10000##b=ion percent in brine
    return b

Ppm_to_Percent(50)
##converting percent to ppm
def Percent_to_Ppm(a):##a=ion percent in brine
    b=a*10000##b=ion concentration in ppm in brine
    return b

Percent_to_Ppm(0.005)

    





#part 1:
vacuum_permeability=1
π=3.14
G=6.674 * 10^-11 
e=2.71828
k=9*10**9
#---------------------
#part2:
def Planks_Fix(y,R,r):
    '''
    This formula is used to calculate the wall pressure of bleeding in the veins
    This function receives the variables y, r, R and returns the variable p
    p:It indicates the pressure of the bleeding wall.
    y:Surface coefficient is surface tension
    R:is the inner radius of the vessel.
    r:is the thickness of the vessel wall.
    '''
    p=2*y*R/r
    return p
#test
zarf=Planks_Fix(10,20,5)
#------------------
def HeatـTransferـCoefficient(k,A,t1,t2,d):
    '''
    

    Parameters
    ----------
    This function receives the variables k, A, t1,t2,d and returns the variable Q
    Q : int
        It indicates the rate of heat transfer.
    k : int
        Heat transfer coefficient
    A : int
        is the heat transfer coefficient.
    t1 : int
        The area of ​​the surface from which heat is transferred.
    t2 : int
        Initial temperature
    d : int
        Secondary temperature

    '''
    Q=k*A*(t2-t1/d)
    return Q
#test
zarf2=HeatـTransferـCoefficient(10,20,11,12,2)
def Power_Factor(i,v):
    '''
    Parameters
    ----------
    This function receives the variables i, v and returns the variable P
    i : int
        It is an electric current that passes through a circuit
    v : int
        It is an electric current that passes through a circuit
    p : int
        It indicates the power used or produced.
    
    '''
    p=i*v
    return p
# test
zarf3=Power_Factor(10,20)
#-----------------
def Electrical_Resistance(v,i):
    '''
    Parameters
    ----------
    This function receives the variables i, v and returns the variable R
    v : int
        Voltage or potential difference between two points
    i : int
        I is the electric current that passes through a circuit
    R : It represents electrical resistance.

    Returns
    -------
    '''
    R =v/i
    return R
# test
zarf4=Electrical_Resistance(10,2)
#-----------------
#part3(convert):
def km_to_hm(km):
    hm = km*10
    return hm
# This function is usef for convert kilometr to hectometr
def hm_to_km(hm):
    km =hm*0.1
    return km
#This function is usef for convert hectometr to kilometr
#-------------------
def km_to_deka(km):
    deka=km*10
    return deka
#This function is usef for convert kilometr to dekametr
def deka_to_km(deka):
    km=deka*0.01
    return km
#This function is usef for convert dekametr to kilometr



#PARET 1
#Define 5 constanr Numbers

# I am polymer and materials engineering so I use constant numbers related to these fields

#part 1

def Constant_Numbers(constant_num):
    
    constants={'Kuhn_Length':1e-9,'e':2.7182818284,'Avogadro':6.02214076e23,'Plank':6.62607015e-34,'K_Boltzman':1.380649e-23,'Chi_parameter':0.5}
    return constants.get(constant_num, None)

print(Constant_Numbers(constant_num='Avogadro'))# we printed them in this way and if it does not define in our function it would return None


#====================================================================================================================================================

#Part 2

#1 related to additive manufacturing

def Print_Time(Print_speed,Volume,Printer_efficiency):
    
    ''' volume(float)= usually Cm3
    
    print_speed(float)= and is speed of print usually cm3 per hour)
    
    Printer_efficience(float)= and between 0 -1)
    
    Retuen:
        
        Print_time=(float)
        the time of print in hour 
    
    '''
    Print_time=float(Volume)/(float(Print_speed)*float(Printer_efficiency))
    return Print_time

n1=float(input('Enter volume'))
n2=float(input('Enter print speed'))
n3=float(input('Enter printer efficiency'))
print(Print_Time(n2, n1, n3))

#------------------------------------------------------------------------------------------------------------------------------------------------



#2
def Lorentz_Lorenz_Constant(n,ro):
    
    ''' n(float)=  refrective index
    ro(float)=density g/cm3 or kg/m3
    
    
    Retuen:
        
        R=(float)
        constant related to light properties of materials(polymers, plstics,...)
    
    '''
    R=float((n**2-1)/(n**2+2)/ro)
    return R

n1=float(input('Enter refrective index of your material'))
n2=float(input('Enter density of your material'))
print(Lorentz_Lorenz_Constant(n1, n2))

#---------------------------------------------------------------------------------------------------------------------------------------------------

#3
import math

def Polymer_Life_Time(E,T,t=1):
    
    ''' E(float)=  activation Energy   J
   T(float)= absolute temperature in K
   t(float)= constant time and usually we consider it 1
    we calculate t and E based on experiences
    
    
    Retuen:
        
        Life_Length=(float)
        the length the life of polymer
    
    '''
    E_kj=E*1000 #because I have math range error I consider this limitation
    exp=E_kj/(Constant_Numbers(constant_num='K_Boltzman')*T)
    if exp>709:
        return float('inf') #if exponent is too large it returns infinity
    
    Life_Length=math.exp(exp)
    return Life_Length

n1=float(input('Enter Activation Energy'))
n2=float(input('Enter Temperature'))
print(Polymer_Life_Time(n1, n2))


#==================================================================================================================================================


#Part 3

#1
def Temperature_C_To_F(temp):
    '''temp=(float) the value of temperature
     
     
     Return:
         float: Temperature in C
     '''
    
    
    tem=(temp*9/5)+32
    return tem


def Temperature_F_To_C(tem):
    '''tem=(float) the value of temperature in C 
     
     
     Return:
         float: Temperature in F
     in this function we convert C to F '''
    
    
    temp=(tem-32)*5/9
    return temp
    
    
   
    
temp=float(input('Enter temerature'))
temp_fahrenheit=Temperature_C_To_F(temp)
print('T in C =',temp_fahrenheit)
temp_Centigrad=Temperature_F_To_C(temp_fahrenheit)
print('T in F=',temp_Centigrad)

#------------------------------------------------------------------------------------------------------------------------------------------------------

#2

def Convert_Viscosity_To_Poise(pa_s):
    '''Pa_s=(float) the viscosity in Pa.S
     
     
     Return:
         float: Viscosity in poise
     '''
    
    
    poise=pa_s*10
    return poise


def Convert_Viscosity_To_Pas(poise):
    '''poise=(float) the viscosity in poise
     
     
     Return:
         float: Viscosity in pas
     '''
    
    
    pa_s=poise/10
    return pa_s
    
    
    
v=float(input('Enter viscosity'))
v1=Convert_Viscosity_To_Poise(v)
print('viscosity  in poise=',v1)
v2=Convert_Viscosity_To_Pas(v1)
print('viscosity  in Pa.s=',v2)





# A1_atena_ jalilzadeh
#َArchitecture







#1(5 ta sabet ....)
#__________________________________________________________________________________________


#Building importance factor
I1=1.4
I2=1.2
I3=1
I4=0.8
#Length to width ratio of works
gold=618.1



#2(3 ta tabe...)
#__________________________________________________________________________________________

 
#meter1

def Meter(tedad,tol,arz,ertefa):
    '''
    

    Parameters
    ----------
    tedad : int
       Similar number of what we want
    tol : int
        the length
    arz : int
        width
    ertefa : int
        Height

    Returns=area
    -------
    TYPE=int 
       

    '''
    
    joz1=(tedad*tol*arz*ertefa)
    return joz1
    joz2=(tedad*tol*arz*ertefa)
    return joz2
    joz3=(tedad*tol*arz*ertefa)
    return joz3
    kol=(joz1+joz2+joz3)
    return kol
    
     
meter(2,11, 3, 2)
meter(1, 2, 5, 10)
meter(10, 5, 5, 4)   




#Determining the weight of one meter of rebar2

def Weight(R):
    '''
    

    Parameters
    ----------
    R : int
        Rebar diameter

    Returns
    -------
    w : int
        The weight of one meter of rebar with a given diameter

    ''' 
    
    w=R**2/162
    return w

Weight(18)




# rebar length3

def Rebar(kham,kaver,dahane):
    '''
    
 
    Parameters
    ----------
    kham : 
        A part of the rebar that bends for connection
    kaver : 
        Two parts of rebar that are placed together
    dahane : 
        The distance between two columns

    Returns
    -------
    l : 
        Calculate the required length of rebar

    '''
    l=((kham*2)+(kaver*2)-dahane)
    return l

rebar(2, 3, 6)


#____________________________________________________________________

#3(2 ta tabe converto)


def Cm_to_m(cm):
    m=cm/100
    return m
    

cm_to_m(100)

 

def M_to_cm(m):
    cm=m*100
    return cm



def M_t0_in:
    In=m*0.0254
    return In


def In_to_m :
    m=In/0.0254
    return m
    



'''
part 1: Constant numbers

'''
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
'''
Planck's constant with symbol
h is a physical constant and the main characteristic of the mathematical formulas of quantum mechanics.
This constant describes the behavior of particles and waves at the atomic scale and the particle aspect of light.
'''
#______________________________________________________________________________________________________
'''
part 2: Functions
'''
#F1________________________________________________________________________________________________
def Fick_Sec_Thin(Thickness,Diffusion_coefficient,Time,Thin_layer_Consistency,Position,Thin_Layer_Metal,Second_metal):
    '''
    Fick's second law predicts how diffusion causes the concentration to change with respect to time.
    In the case where a thin film of a different metal is placed between two thick films of another metal,
    and due to heat and a certain time, the amount of penetration of the first metal into the second metal will be investigated.
    Parameters
    ----------
    Thickness : float 
        The thickness of the thin layer metal is used. use cm for input
    Diffusion_coefficient : float
        It shows the diffusion coefficient of the metal placed in the middle compared to the metal on its side. use cm^2/seconds for input
    Time : float
        It shows the duration of the diffusion check. use seconds for input
    Thin_layer_Consistency : float
        It shows the initial concentration of thin metal. use gr/cm^3 for input
    Position :  float
        Indicates the location required to check the concentration. use cm for input
    Thin_Layer_Metal : str
        What is the material of the thin layer used?
    Second_metal : str
        What is the metal material in which the primary metal should penetrate?.

    Returns
    -------
    C_x_t : float
        It shows the concentration of the metal separated from the thin layer and moved in a certain time in the desired location. output will be gr/cm^3
    '''
    import math
    pi=3.14
    if Thin_Layer_Metal=='Cu'and Second_metal=='Ni':
        Diffusion_coefficient=0.2698 # @Temprature=1000 K
    if Thin_Layer_Metal=='Cr'and Second_metal=='Ni':
        Diffusion_coefficient=0.0299 # @Temprature=1000 K
    C_x_t=((Thickness*Thin_layer_Consistency)/(2*(pi*Diffusion_coefficient*Time)**(0.5)))*math.exp((-(Position)**2)/(4*Diffusion_coefficient*Time))
    return C_x_t
'''
in tabe kheili ja kar dare bara khode Diffusion Coefficient mishe ye tabe tarif kard o edame dad vali goftm dar hadi k data peida kardm flan erae bdm ta badan sare vqt kamelesh knm
'''
#F2________________________________________________________________________________________________
def Final_Temp_Irreversible_Adiabatic(Initial_temperature,External_pressure,Internal_pressure,C_V,C_P,R,Unit_of_measurement,Number_of_gas_atoms):
    '''
    Parameters
    ----------
  Initial_temperature : float
      Initial temperature of an ideal gas.
  External_pressure : float
      The pressure that enters the system from the environment.
  Internal_pressure : float
      The pressure that enters the system wall from the inside.
  C_V : float
      The ideal gas constant-pressure specific heat.
  C_P : float
      The ideal gas constant-pressure specific heat..
  R : float
      The molar gas constant or ideal gas constant.
  Unit_of_measurement : str
      DESCRIPTION.
  Number_of_gas_atoms : str
      DESCRIPTION.

  Returns
  -------
  Final_Temp : float
      Outputs the final temperature if the system operates adiabatically irreversible.
    '''
    if Unit_of_measurement=='Jouls':
        R=8.314
        if Number_of_gas_atoms=='one':
            C_V=1.5*R
            C_P=2.5*R
        if Number_of_gas_atoms=='two':
            C_V=2.5*R
            C_P=3.5*R
        if Number_of_gas_atoms=='multi':
            C_V=3.5*R
            C_P=4.5*R
    elif Unit_of_measurement=='Calories':
        R=1.98
        if Number_of_gas_atoms=='one':
            C_V=1.5*R
            C_P=2.5*R
        if Number_of_gas_atoms=='two':
            C_V=2.5*R
            C_P=3.5*R
        if Number_of_gas_atoms=='multi':
            C_V=3.5*R
            C_P=4.5*R   
    Final_Temp=Initial_temperature*((C_V+(External_pressure/Internal_pressure)*R)/C_P)
    return Final_Temp
#F3________________________________________________________________________________________________
def Atomic_Percentage(n_1,n_2,n_3,n_4,Entry,output,m_1,m_2,m_3,m_4,M_1,M_2,M_3,M_4,w_1,w_2,w_3,w_4):
    '''
    Parameters
    ----------
    n_1 : float
        The number of moles of the first atom.
    n_2 : float
        The number of moles of the second atom.
    n_3 : float
        The number of moles of the third atom.
    n_4 : float
        The number of moles of the fourth atom.
    Entry : str
        what type of entery you have?
    output : str
        The atomic percentage is based on which atom?
    m_1 : float
        The mass of the first atom.
    m_2 : float
        The mass of the second atom.
    m_3 : float
        The mass of the third atom.
    m_4 : float
        The mass of the fourth atom.
    M_1 : float
        The atomic mass of the first atom.
    M_2 : float
        The atomic mass of the second atom.
    M_3 : float
        The atomic mass of the third atom.
    M_4 : float
        The atomic mass of the fourth atom.
    w_1 : float
        The weight percentage of the first atom.
    w_2 : float
        The weight percentage of the second atom.
    w_3 : float
        The weight percentage of the third atom.
    w_4 : float
        The weight percentage of the fourth atom.

    Returns
    -------
    float
        It will return the atomic percentage based on the given inputs.
    '''
    if Entry=='mole':
        if output=='AP_1': #atomic percent for first atom
            AP_1=(n_1)/(n_1+n_2+n_3+n_4)
            return AP_1
        elif output=='AP_2': #atomic percent for second atom
            AP_2=(n_2)/(n_1+n_2+n_3+n_4)
            return AP_2
        elif output=='AP_3': #atomic percent for third atom
            AP_3=(n_3)/(n_1+n_2+n_3+n_4)
            return AP_3
        elif output=='Ap_4': #atomic percent for fourth atom
            AP_4=(n_4)/(n_1+n_2+n_3+n_4)
            return AP_4
    if Entry=='mass':
        if output=='AP_1':
            AP_1=(m_1/M_1)/((m_1/M_1)+(m_2/M_2)+(m_3/M_3)+(m_4/M_4))
            return AP_1
        if output=='AP_2':
            AP_2=(m_2/M_2)/((m_1/M_1)+(m_2/M_2)+(m_3/M_3)+(m_4/M_4))
            return AP_2
        if output=='AP_3':
            AP_3=(m_3/M_3)/((m_1/M_1)+(m_2/M_2)+(m_3/M_3)+(m_4/M_4))
            return AP_3
        if output=='AP_4':
            AP_4=(m_4/M_4)/((m_1/M_1)+(m_2/M_2)+(m_3/M_3)+(m_4/M_4))
            return AP_4
    if Entry=='weight':
        if output=='AP_1':
            AP_1=(w_1/M_1)/((w_1/M_1)+(w_2/M_2)+(w_3/M_3)+(w_4/M_4))
            return AP_1
        if output=='AP_2':
            AP_2=(w_2/M_2)/((w_1/M_1)+(w_2/M_2)+(w_3/M_3)+(w_4/M_4))
            return AP_2
        if output=='AP_3':
            AP_3=(w_3/M_3)/((w_1/M_1)+(w_2/M_2)+(w_3/M_3)+(w_4/M_4))
            return AP_3
        if output=='AP_4':
            AP_4=(w_4/M_4)/((w_1/M_1)+(w_2/M_2)+(w_3/M_3)+(w_4/M_4))
            return AP_4
'''
part 3: Convertors
'''
#C1________________________________________________________________________________________________
def Hour_to_Second(Hour):
    Second=Hour*3600
    return Second

def Second_to_Hour(Second):
    Hour=Second/3600
    return Hour
#C2________________________________________________________________________________________________
def Joules_to_Calories(Joules):
    Calories=Joules/4.2
    return Calories

def Calories_to_Joules(Calories):
    Joules=Calories*4.2
    return Joules




#-----------------------------------------------------------------------------

import math

#------------------1st--------------------------------------------------------

R = 8.314                        # Gas constant        --> j/mol*k

a = 6.0232 * (10**23)            # Avogadro constant   --> 1/mol

k = 1.38054 * (10**23)           # Boltzmann constant  --> joules/degree

F = 2306                         # Faraday constant    --> calories/volt*mol

h = 6626068 * (10**34)           # Planck constant     --> J*s

c = 4.18                         # Specific heat capacity of liquid water 

m = 9.1093837015 * (10**(-34))   # Electron mass --> kg

#------------------2nd--------------------------------------------------------

# 1st function calculates Latice parameter of crystal structure.

# FCC(face centered cubic) 
# BCC(body centered cubic) 
# SC(simple cubic) 
# HCP(hexagonal close pack)
# DC(diamond cubic)

def Latice_Parameter(structure, r):
    
    '''
    this function calculates Latice parameter based on crystal structure.
    
    Parametrs
    ----------
    1. structure : str
    structure includes crystal structures such as FCC, BCC, SC, HCP, and DC.
    
    FCC: face centered cubic
    BCC: body centered cubic
    SC:  simple cubic
    HCP: hexagonal close pack
    DC:  diamond cubic
    
    2. r: float ---> (nm)
    r represents atomic radius.

    Returns --> latice parameter (a)
        
    '''
    
    
    
    if structure == 'FCC' or  'face centered cubic':
        a = (4*r)/math.sqrt(2)
        
    elif structure == 'BCC' or  'body centered cubic':
        a = (4*r)/math.sqrt(3)
     
    elif structure == 'SC' or  'simple cubic':
        a = 2*r  
     
    elif structure == 'HCP' or  'hexagonal close pack':
        a = 2*r  
                            
    elif structure == 'DC 'or  'diamond cubic':
        a = (8*r)/math.sqrt(3)                        
    
    return a  


# call function----------------------------------------------------------------

Latice_Parameter('FCC', 0.1431 )  # Aluminium
Latice_Parameter('face centered cubic', 0.1387 )  # Platinum


Latice_Parameter('BCC', 0.1249 )  # Chromium
Latice_Parameter('body centered cubic', 0.1371 )  # Tungsten


Latice_Parameter('HCP', 0.1332 )  # Zinc
Latice_Parameter('hexagonal close pack', 0.1445 )  # Titanium


parameter = Latice_Parameter('FCC', 0.1750 )  # Lead
print(parameter)

#------------------------------------------------------------------------------   

# 2nd function calculates fracture toughness 
# based on applied stress and crack length and location

def Fracture_Toughness(s, c, location):
    
    '''
    this function calculates fracture toughness 
    based on applied stress, crack length, and location
    
     fracture toghness formula:      
         K1C =  a*s* (math.sqrt(3.14* c))
         
    where:

        K1C is the fracture toughness,
        a is a geometric factor,
        s is the applied stress,
        c is the crack length.

         
    Parameters
    ----------- 
    1. s: int ---> (MPa)
       applied stress 
       
    2. c: float ---> (m) 
       crack length
       
    3. location: str
       represents the location of crack

    Returns --> K1C (fracture toughness)
    
    '''
    
    if location == 'surface' or 'side':
        a = 1.12
        K1C =  a*s* (math.sqrt(3.14* c))
        
    elif location == 'centeral':
        a = 1
        c = c/2
        K1C =  a*s* (math.sqrt(3.14* c))
        
        
    return K1C


# call function----------------------------------------------------------------

Fracture_Toughness(10, 0.0796, 'surface')

#------------------------------------------------------------------------------       

# 3rd function calculates wear rate.

def Wear_Rate(v, f, s):
    
    '''
    this function calculates abrasive wear rate
    
    Parametrs
    ---------
    1. v : float ---> (mm3)
    v represents wear volume
    
    2. f: int ---> (nm) ---> (N)
    f represents normal load.
    
    3. s: float ---> (m)
    s represents sliding distance.
    
    Returns --> wear rate (w)
        
    '''
    
    w = v/(f*s)
    
    return w

# call function----------------------------------------------------------------

Wear_Rate(0.0006, 20, 100)

#------------------------------------------------------------------------------

# 4th function calculates Vickers hardness.

def Vicker_Hardness_Calculation (d1,d2,p): 
    
    '''
    this function is utilized for calculating Vickers hardness.
    in this method, pyramid shape indentor is used.
    
    Parameters
    -----------
    1. d1: float
    d1 represents diameter of impress of indentor.
    
    2. d2 : float
    d2 represents diameter of impress of indentor
    
    3. p : int
    applied load 
    
    Returns --> VHN (Vickers Hardness Number)
    '''
    
    d= (d1+d2)/2
    VHN = (1.854*p)/(d**2)

    return VHN

#------------------------------------------------------------------------------

## 5th function calculates Brinell hardness.

def Brinell_Hardness_Calculation (d1,d2,D,p): 
    
    '''
    this function is utilized for calculating Brinell hardness.
    characterizes the indentation hardness of materials through the
    scale of penetration of an indenter, loaded on a material test-piece. 
    It is one of several definitions of hardness in materials science.
    In this method, steel ball indentor is used.
    
     Parameters
     ----------- 
            1. d1: float
            d1 represents diameter of impress of indentor.
            
            2. d2 : float
            d2 represents diameter of impress of indentor
            
            3. D : float
            diameter of indenter
            
            3. p : int
            applied load 
     
            Returns --> BHN (Brinell Hardness Number)
    '''
    
    d = (d1+d2)/2
    BHN = (2*p)/(3.14*D)(D-(math.sqrt(D**2 - d**2)))

    return BHN

#------------------3rd---------------------------------------------------------

# 1st function converts torr to pascal and vice versa.

def Torr_to_Pascal(torr):
    
    '''
    this function converts torr to pascal
    '''
    
    
    pa = (torr*760)/101325
    
    return pa
#------------------------------------------------------------------------------

def Pascal_to_Torr(pa):
    
    '''
    this function converts pascal to torr
    '''
    
    torr = (pa*101325)/760
    
    return torr
#------------------------------------------------------------------------------

# 2nd function converts joules to calories and vice versa.

def Joules_to_Calories(j):
    
    '''
    this function converts joules to calories
    '''
    
    
    cal = j / 4.184
    
    return cal
#------------------------------------------------------------------------------

def Calories_to_Joules(cal):
    
    '''
    this function converts calories to joules
    '''
    
    j = 0.239006 * cal
    
    return j
#------------------------------------------------------------------------------

# 3rd function converts atmosphere to pascal and vice versa.

def Atmosphere_to_Pascal(atm):
    
    '''
    this function converts atmosphere to pascal
    '''
    
    
    pa = atm * 101325
    
    return pa
#------------------------------------------------------------------------------

def Pascal_to_Atmosphere(pa):
    
    '''
    this function converts pascal to atmosphere
    '''
    
    atm = pa/101325
    
    return atm


#------------------------------------------------------------------------------

# 4th function converts miller index to miller_brove index and vice versa

# In crystallography, crystal structure is a description of the ordered 
# arrangement of atoms, ions, or molecules in a crystalline material.
# unite cell repeat along the principal directions of three-dimensional space.
# but there are 2 patterns in HCP structure. 
# one of them includes 3 component which called miller index.
# another involves 4 component which called miller_brove index. 
# mentioned indexes determine crystal directions.


def Miller_to_Millerbrove(u,v,w):
    
    '''
       this function converts miller index to miller_brove index

     parameters: (miller indexes)
     ---------------------------  
        1. u: int
        Intersection with axis a1
        
        2. v: int
        Intersection with axis a2
        
        3. w: int
        Intersection with axis z
        
    Returns --> (miller_brove indexes)
            
       1. l: int
       Intersection with axis a1
       
       2. m: int
       Intersection with axis a2
       
       3. n: int
       Intersection with axis a3
       
       4. o: int
       Intersection with axis z
      
  '''
  
    l = ((2*u)-v)/3
    m = ((2*v)-u)/3
    n = (-1)*(l+m)
    o = w
    
    return l,m,n,o
  
#------------------------------------------------------------------------------

def Millerbrove_to_Miller(l,m,n,o):

    '''
       this function converts miller_brove index to miller index

    Parameters: (miller_brove indexes)
    -----------------------------------
        1. l: int
       Intersection with axis a1
       
       2. m: int
       Intersection with axis a2
       
       3. n: int
       Intersection with axis a3
       
       4. o: int
       Intersection with axis z
       
      Returns --> (miller indexes)
             
        1. u: int
        Intersection with axis a1
        
        2. v: int
        Intersection with axis a2
        
        3. w: int
        Intersection with axis z
  '''
    u = (2*m) + l
    v = (2*l) + m
    w = o
    
    
    return u,v,w






# -*- coding: utf-8 -*-
"""
A1_Zahra_Afshar

send to ai.2024.pilehvar@gmail.com

@author: Zahra Afshar

"""

#1:
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








#2:
    #2-1:
def Gibs_free_energy(H0,T,S0):
    '''
    Parameters
    ----------
    H0 : float
        Enthalpy of material in the pressure of 1atm in pure state. The number should be in terms of joul or kilojoul. It must have the same unit as entropy.
    T : float
        Temperature of material. The temperature should be in terms of Kelvin.
    S0 : float
        Entropy of material. The number should be in terms of joul or kilojoul. It must have the same unit as enthalpy.

    Returns
    -------
    G : float
        This function give us the amount of Gibs free energy of the material at the pressure of 1atm. The terms of this function is joul or kilojoul.

    '''
    G0=H0-T*S0
    return G0
Gibs_free_energy(25,10,2)




    #2-2:
def Hardness_vickers(F,d):
    '''
    Parameters
    ----------
    F : float
        Applied force in terms of kilogram.
    d : float
        Medium diameter of indentor effect in terms of milimeter.

    Returns
    -------
    HV : Vickers hardness
        This function give us the amount of hardness of material that evaluated by vickers hardness tester. The terms of this is Kg/mm2.

    '''
    HV=1.854*F/(d**2)
    return HV
Hardness_vickers(20, 2)






    #2-3:
def Wear_rate(V,F,S):
    '''
    Parameters
    ----------
    V : float
        lost volume during the wear test in terms of cubic milimeters (mm3)
    F : float
        Applied force in terms of newton.
    S : float
        Sliding distance in terms of meter.

    Returns
    -------
    K : Wear rate
        This function is for claculation of wear rate if material after wear test. It is in terms of (10-6 mm3/N.m) 
        

    '''
    K=V/(F*S)
    return K
Wear_rate(120, 20, 100)




   #2-4:
def Lattice_Parameter(r,structure):
    '''
    Parameters
    ----------
    r : float
        r is atomic radius of material. It is in terms of angestrom.
    structure : str
        Structure of material is face center cubic or body center cubic.

    Returns
    -------
    a : float
        a is lattice parameter of material in unit cell. It is in terms of angestrom.

    '''
    if structure=='fcc':
        a=4*r/(2**0.5)
    if structure=='bcc':
        a=4*r/(3**0.5)
    return a

Lattice_Parameter(2, 'bcc')

  







  
#3:
#3-1:
def Torr_Pa(Torr):
    '''
    This function convert the amount of pressure in terms of Torr to Paskal.    
    '''
    Pa=Torr*133.322
    return Pa
Torr_Pa(2)


def Pa_Torr(Pa):
    '''
    This function convert the amount of pressure in terms of paskal to torr.
    '''
    Torr=Pa/133.322
    return Torr
Pa_Torr(266.644)


#3-2:
def kgf_dyn(kgf):
    '''
    This function convert the amount of force in terms of kilogramforce to dyn.
    '''
    dyn=kgf*980665
    return dyn
kgf_dyn(3)


def dyn_kgf(dyn):
    '''
    This function convert the amount of force in terms of dyn to kilogramforce.
    '''
    kgf=dyn/980665
    return kgf
dyn_kgf(2941995)

    
#A1_Project
#Part 1

Boltzmann_Constant=1.380649*(10**(-23)) #joule per kelvin (J/K)
Avogadro_Constant=6.02214*(10**23) #per moles (mol-1)
Faraday_Constant=96485.3399 #coulombs per mole of electrons (C/mol)
Planck_Constant=6.62607015*(10**(-34)) #joule second (J.s)
Elementary_Charge=1.602176634*(10**(-19)) #coulombs (C)
Light_Speed_Constant=299792458 #meters per second (m/s)

#------------------------------
#Part 2

def Thermal_Voltage(T1,T2):
    
    """
    
    Parameters
    ----------
    T1 : float
        Temperature 1 (in Kelvin).
    T2 : float
        Temperature 2 (in Kelvin).

    Returns
    -------
    V_T : float
        Thermal Voltage is voltages created by the junction of dissimilar metals when a temperature difference exists between these junctions.

    """
    
    V_T=(Boltzmann_Constant*(abs(T2-T1))/Elementary_Charge)
    return V_T
        

def Number_of_Particles(n):
    
    """
    
    Parameters
    ----------
    n : float
        Number of moles.

    Returns
    -------
    N : float
        This function is used to calculate the number of particles or atoms in a substance.

    """
    
    N=n*Avogadro_Constant
    return N

def Faraday_Law_of_Electrolysis(Q,M,z):
    """
    
    Parameters
    ----------
    Q : float
        Total charge (in coulombs).
    M : float
        Molar mass of the substance (in grams per mole).
    z : float
        Number of moles of electrons transferred per mole of substance.

    Returns
    -------
    m : float
        Mass of substance deposited or liberated (in grams)

    """
    
    m=(Q*M)/(z*Faraday_Constant)
    return m
    
def Photon_Energy(f):
    """

    Parameters
    ----------
    f : float
        Frequency of the photon (in hertz).

    Returns
    -------
    E : float
        Energy of the photon (in Joule).

    """
    
    E=Planck_Constant*f
    return E
    
#------------------------------
#Part 3

def Calories_to_Joules(cal):
    """

    Parameters
    ----------
    cal : float
        Calories.

    Returns
    -------
    J : float
        Converts calories to joules.

    """
    
    J=4.184*cal
    return J

def Joules_to_Calories(J):
    """

    Parameters
    ----------
    J : float
        Joules.

    Returns
    -------
    cal : float
        Converts joules to calories.

    """
    
    cal=J/4.184
    return cal

def Bar_to_Pascal(bar):
    """

    Parameters
    ----------
    bar : float
        bar.

    Returns
    -------
    Pa : float
        Converts bar to pascal.

    """
    
    Pa=bar*(10**(-5))
    return Pa

def Pascal_to_Bar(Pa):
    """

    Parameters
    ----------
    Pa : float
        Pascal.

    Returns
    -------
    bar : float
        Converts pascal to bar.

    """
    
    bar=Pa*(10**5)
    return bar



    


















'''

#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================
#========================

#========================
#A2

'''    


def Diabetes_Dataset(f, work):
    """
    Reads a diabetes dataset from a CSV file

    Parameters:
        f (str): The file path of the CSV file containing the diabetes dataset.
        work (str): The task to perform. Supported values are:
            - 'has_diabetes': Counts the number of individuals who have diabetes.
            - 'percent_has_diabetes_25': Calculates the percentage of individuals under 25 who have diabetes.
            - 'percent_has_diabetes_25_to_30': Calculates the percentage of individuals between 25 to 30 who have diabetes.
            - 'percent_has_diabetes_30_to_40': Calculates the percentage of individuals between 30 to 40 who have diabetes.
            - 'percent_has_diabetes_40_and_50': Calculates the percentage of individuals between 40 to 50 who have diabetes.
            - 'percent_has_diabetes_20_to_40': Calculates the percentage of individuals between 20 to 40 who have diabetes.
            - 'percent_has_diabetes_30_to_50': Calculates the percentage of individuals between 30 to 50 who have diabetes.
            - 'percent_has_diabetes_50_80': Calculates the percentage of individuals between 50 to 80 who have diabetes.
            - 'rel_bmi_to_diabetes_30_to_40': Calculates the percentage of individuals with BMI between 30 to 40 who have diabetes.
            - 'rel_bmi_to_diabetes_20_to_30': Calculates the percentage of individuals with BMI between 20 to 30 who have diabetes.
            - 'histo': Plots a histogram of ages.
            - 'max_age': Finds the maximum age in the dataset.
            - 'min_age': Finds the minimum age in the dataset.
            - 'max_age_has_diabetes': Finds the individuals with the maximum age who have diabetes.
            - 'min_age_has_diabetes': Finds the individuals with the minimum age who have diabetes.
    """
    read_data = pd.read_csv(f)

    if work == 'has_diabetes':
        all_outcome = read_data['Outcome']
        has_diabetes = []
        for i in all_outcome:
            if i == 1:
                has_diabetes.append(i)
        text = f'all Outcome is {len(all_outcome)} and has diabetes is {len(has_diabetes)}'
        return text
    elif work == 'percent_has_diabetes_25':
        has_diabetes_25 = len(read_data[(read_data['Age'] <= 25) & (read_data['Outcome'] == 1)])
        less_then_25 = len(
            read_data[(read_data['Age'] <= 25) & ((read_data['Outcome'] == 1) | (read_data['Outcome'] == 0))])
        percent = has_diabetes_25 * (100 / less_then_25)
        return percent
    elif work == 'percent_has_diabetes_25_to_30':
        between_25_to_30 = read_data[((read_data['Age'] > 25) & (read_data['Age'] <= 30)) & (read_data['Outcome'] == 1)]
        all_outcome_25_to_30 = read_data[((read_data['Age'] > 25) & (read_data['Age'] <= 30)) &
                                         ((read_data['Outcome'] == 1) | read_data['Outcome'] == 0)]
        percent = len(between_25_to_30) * (100 / len(all_outcome_25_to_30))
    elif work == 'percent_has_diabetes_30_to_40':
        between_30_to_40 = read_data[((read_data['Age'] > 30) & (read_data['Age'] <= 40)) & (read_data['Outcome'] == 1)]
        all_outcome_30_to_40 = read_data[((read_data['Age'] > 30) & (read_data['Age'] <= 40)) &
                                         ((read_data['Outcome'] == 1) | (read_data['Outcome'] == 0))]
        percent = len(between_30_to_40) * (100 / len(all_outcome_30_to_40))
        return percent
    elif work == 'percent_has_diabetes_40_to_50':
        has_diabetes_40_to_50 = read_data[((read_data['Age'] > 40) & (read_data['Age'] <= 50)) &
                                          (read_data['Outcome'] == 1)]
        all_outcome_40_to_50 = read_data[((read_data['Age'] > 40) & (read_data['Age'] <= 50)) &
                                         ((read_data['Outcome'] == 1) | (read_data['Outcome'] == 0))]
        percent = len(has_diabetes_40_to_50) * (100 / len(all_outcome_40_to_50))
        return percent

    elif work == 'percent_has_diabetes_20_to_40':
        age_groups = [(18, 30)]
        for age_group in age_groups:
            min_age, max_age = age_group
            data = read_data[(read_data['Age'] >= min_age) & (read_data['Age'] <= max_age)]
            has_diabetes = len(data[data['Outcome'] == 1])
            all_data = len(data)
            percent = has_diabetes * (100 / all_data)
            return percent
    elif work == 'percent_has_diabetes_30_to_50':
        age_groups = [(30, 50)]
        for age_group in age_groups:
            min_age, max_age = age_group
            data = read_data[(read_data['Age'] >= min_age) & (read_data['Age'] <= max_age)]
            has_diabetes = data[data['Outcome'] == 1]
            percent = len(has_diabetes) * (100 / len(data))
            return percent
    elif work == 'percent_has_diabetes_50_80':
        age_groups = [(50, 80)]
        for age_group in age_groups:
            min_age, max_age = age_group
            data = read_data[(read_data['Age'] >= min_age) & (read_data['Age'] <= max_age)]
            has_diabetes = data[data['Outcome'] == 1]
            percent = len(has_diabetes) * (100 / len(data))
            return percent
    elif work == 'rel_bmi_to_diabetes_30_to_40':
        bmi_groups = [(30, 40)]
        for bmi_group in bmi_groups:
            min_bmi, max_bmi = bmi_group
            data = read_data[(read_data['BMI'] >= min_bmi) & (read_data['BMI'] <= max_bmi)]
            has_diabetes = data[data['Outcome'] == 1]
            percent = len(has_diabetes) * (100 / len(data))
            return percent
    elif work == 'rel_bmi_to_diabetes_20_to_30':
        bmi_groups = [(20, 30)]
        for bmi_group in bmi_groups:
            min_bmi, max_bmi = bmi_group
            data = read_data[(read_data['BMI'] >= min_bmi) & (read_data['BMI'] < max_bmi)]
            has_diabetes = data[data['Outcome'] == 1]
            percent = len(has_diabetes) * (100 / len(data))
            return percent
    elif work == 'histo':
        ages = read_data['Age']
        plt.hist(ages, bins=20, color='skyblue', edgecolor='green')
        plt.show()
    elif work == 'max_age':
        age = read_data['Age']
        max_age = np.max(age)
        return max_age
    elif work == 'min_age':
        age = read_data['Age']
        min_age = np.min(age)
        return min_age
    elif work == 'max_age_has_diabetes':
        max_ages = read_data['Age'].max()
        has_diabetes = read_data[((read_data['Age'] == max_ages) & (read_data['Outcome'] == 1))]
        return has_diabetes
    elif work == 'min_age_has_diabetes':
        min_ages = read_data['Age'].min()
        has_diabetes = read_data[((read_data['Age'] == min_ages) & (read_data['Outcome'] == 1))]
        return has_diabetes
    else:
        raise Exception('invalid command')


def Income_Developer(file, work):
    """
    Reads income data from a CSV file and displays graphs based on the specified task.

    Parameters:
        file (str): The file path of the CSV file containing the income dataset.
        work (str): The task to perform. Supported values are:

        - 'plot_data': Displays a line plot of income over age for Python and JavaScript developers.
        - 'bar_data': Displays a bar plot of income over age for Python and JavaScript developers.
        - 'max_salary_data': Displays a bar plot showing the maximum income for Python and JavaScript developers.
        - 'plot_bar_data': Displays both bar and line plots of income over age for Python and JavaScript developers.
        - 'max_salary_data_by_age': Displays a bar plot showing the maximum income for Python and JavaScript developers based on age.
        - 'alph_data': Displays a bar plot with different alpha values for Python and JavaScript developers.
        - 'show_by_side_by_side_data': Displays two bar plots side by side for Python and JavaScript developers with age on the x-axis and income on the y-axis.

    """

    # read file
    data = pd.read_csv(file)

    # style plt
    plt.style.use('fast')

    # get data
    python_data = data[(data['Language'] == 'python')]
    js_data = data[(data['Language'] == 'js')]

    # get age data
    python_ages = python_data['Age']
    js_ages = js_data['Age']
    all_age = data['Age']
    ages_x = np.arange(18, 49)

    # get income data
    python_income = python_data['Income']
    js_income = js_data['Income']

    # get max income data
    max_income_python = python_data['Income'].max()
    max_income_js = js_data['Income'].max()

    # convert data to array
    np_python_ages = np.array([python_ages])
    np_js_ages = np.array([js_ages])

    if work == 'plot_data':
        # show plot data
        plt.plot(python_ages, python_income, color='skyblue', label='Python')
        plt.plot(js_ages, js_income, color='green', label='js')
        plt.legend()
        plt.show()
    elif work == 'bar_data':
        # show bar data
        plt.bar(python_ages, python_income, color='skyblue', label='Python')
        plt.bar(js_ages, js_income, color='green', label='js')
        plt.legend()
        plt.show()
    elif work == 'max_salary_data':
        plt.bar(['python', 'js'], [max_income_python, max_income_js], color=['skyblue', 'green'],
                label=['python', 'js'])
        plt.legend()
        plt.show()
    elif work == 'plot_bar_data':
        plt.bar(python_ages, python_income, color='skyblue', label='Python')
        plt.bar(js_ages, js_income, color='green', label='js')
        plt.plot(python_ages, python_income, color='blue', label='python')
        plt.plot(js_ages, js_income, color='yellow', label='js')
        plt.legend()
        plt.show()
    elif work == 'max_salary_data_by_age':
        python_row = python_data[python_data['Income'] == max_income_python]
        js_row = js_data[js_data['Income'] == max_income_js]
        age_py = python_row['Age']
        age_js = js_row['Age']
        income_py = python_row['Income']
        income_js = js_row['Income']
        plt.bar(age_py, income_py, color='skyblue', label='Python')
        plt.bar(age_js, income_js, color='green', label='js')
        plt.legend()
        plt.show()
    elif work == 'alph_data':
        plt.bar(python_ages, python_income, color='skyblue', label='python', alpha=0.8)
        plt.bar(js_ages, js_income, color='green', label='js', alpha=0.4)
        plt.legend()
        plt.show()
    elif work == 'show_by_side_by_side_data':
        my_width = 0.6
        age = np.arange(len(np.arange(18, 49)))
        plt.bar(age + my_width, python_income, color='skyblue', label='python', width=0.4)
        plt.bar(age, js_income, color='green', label='js', width=0.4)
        plt.xticks(ticks=age, labels=ages_x)
        plt.xlabel('age')
        plt.ylabel('salary')
        plt.legend()
        plt.show()
    else:
        raise Exception('invalid command')


# f1 = os.path.join('diabetes.csv')
f2 = os.path.join('income_developer.csv')
Income_Developer(f2, 'show_by_side_by_side_data')
# Diabetes_Dataset(f1, '')



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from tensile_test_data import F_numeric_values, DL_numeric_values
import math

tensile_test_data = {'DL': DL_numeric_values, 'F': F_numeric_values}
df = pd.DataFrame(tensile_test_data)


def Strain_Stress(df,operation):
    '''
    
    
    This function gets data and an operation .
    It plots Stress-Strain curve if the oepration is plot 
    and finds the UTS value (which is the ultimate tensile strength) otherwise.
    ------------------------------
    Parameters
    ----------
    df : DataFrame
       It has 2 columns: DL(which is length in mm) & F (which is the force in N).
    operation :
       It tells the function to whether Plot the curve or find the UTS valu. 

    Returns
    -------
    The Stress-Strain curve or the amount of UTS
    
    '''
    L0 = 40
    #L0: initial length of the sample
    D0 = 9
    #D0: initial diameter of the sample
    A0 = math.pi / 4 * (D0 ** 2)
    df['e'] = df['DL'] / L0
    df['S'] = df['F'] / A0
    if operation == 'PLOT':
        plt.scatter(df['e'], df['S'])
        plt.xlabel('e')
        plt.ylabel('S')
        plt.title('S vs e Plot')
        plt.grid(True)
        plt.show()
    elif operation == 'UTS':
        return df['S'].max()
    else:
        print("Please enter proper operation")
        return

def LN_S_E(df, operation):
    '''
    

    This function plots the ELASTIC part of the true stress-strain curve
    and can also calculate the Young's Modulus(E).
    ----------
    df : DataFrame
       It has 2 columns: DL(which is length in mm) & F (which is the force in N).
    operation :
       It tells the function to whether Plot the curve or calculate the Young's Modulus(E).

    Returns
    -------
    The elastic part of the true stress-strain curve or the amount of E

    '''
    L0 = 40
    D0 = 9
    A0 = math.pi / 4 * (D0 ** 2)
    df['e'] = df['DL'] / L0
    df['S'] = df['F'] / A0
    df['eps'] = np.log(1 + df['e'])
    df['sig'] = df['S'] * (1 + df['e'])
    filtered_index = df[(df['eps'] >= 0.04) & (df['eps'] <= 0.08)].index
    "the elastic part of the curve is where the true strain(eps) is from 0.04 to 0.08"

    df['selected_ln_eps'] = np.nan

    df.loc[filtered_index, 'selected_ln_eps'] = df.loc[filtered_index, 'eps']
    df['selected_ln_eps'] = np.where(~df['selected_ln_eps'].isna(), np.log(df['selected_ln_eps']), df['selected_ln_eps'])
    df['selected_ln_sig'] = np.nan

    df.loc[filtered_index, 'selected_ln_sig'] = df.loc[filtered_index, 'sig']
    df['selected_ln_sig'] = np.where(~df['selected_ln_sig'].isna(), np.log(df['selected_ln_sig']), df['selected_ln_sig'])


    if operation == 'PLOT':
        plt.scatter(df['selected_ln_eps'].dropna(), df['selected_ln_sig'].dropna())
        plt.xlabel('ln_eps')
        plt.ylabel('ln_sig')
        plt.title('ln(sig) vs ln(eps) Plot')
        plt.grid(True)
        plt.show()
    elif operation == 'YOUNG_MODULUS':
        X = df['selected_ln_eps'].dropna().values.reshape(-1, 1)  # Independent variable
        y = df['selected_ln_sig'].dropna().values.reshape(-1, 1)  # Dependent variable
        model = LinearRegression()
        model.fit(X, y)
        intercept = model.intercept_[0]
        return math.exp(intercept)
        
    else: 
        print("Please enter proper operation")
        return

Strain_Stress(df,'PLOT')
print("UTS is: " , Strain_Stress(df,'UTS'))

LN_S_E(df,'PLOT')
print("Young's Modulus is: " , LN_S_E(df,'YOUNG_MODULUS'))
def Oxygen_Heat_Capacity_Analysis(file_path):
    '''
    This function reads the temperature and heat capacity information
    from the Excel file and then calculates the enthalpy and entropy values
    and draws their graphs.

    Parameters:
    T: Oxygen temperature
    Cp: Heat capacity at constant pressure in different oxygen temperature ranges
    
    Return:
    The values of enthalpy, entropy and their graph according to temperature
    Heat capacity graph according to temperature
    '''

    data = pd.read_excel(file_path)
    
    print("value of T column")
    print(data["T"])
    
    print("value of cp column")
    print(data["Cp"])
    
    data['Enthalpy'] = data['Cp'].cumsum()
    data['Entropy'] = data['Enthalpy'] / data['T']

    fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    axs[0].plot(data['T'], data['Cp'])
    axs[0].set_xlabel('Temperature (T)')
    axs[0].set_ylabel('Heat Capacity (Cp)')
    axs[0].set_title('Heat Capacity vs Temperature')

    axs[1].plot(data['T'], data['Enthalpy'], label='Enthalpy')
    axs[1].plot(data['T'], data['Entropy'], label='Entropy')
    axs[1].set_xlabel('Temperature (T)')
    axs[1].set_ylabel('Enthalpy and Entropy')
    axs[1].set_title('Enthalpy and Entropy vs Temperature')
    axs[1].legend()

    plt.tight_layout()
    plt.show()

file_path = 'C:/Data/Oxygen_Heat_Capacity.xlsx'  
Oxygen_Heat_Capacity_Analysis(file_path)

def Stress_Strain_Analysis(file_path, D0, L0):
    '''
    This function uses the data file
    that contains length and force, calculates the engineering, true
    and yielding stress and strain and also draws a graph of these.
    
    Parameters:
    D0(mm): First Qatar to calculate stress
    L0(mm): First Length to canculate strain
    F(N): The force applied to the object during this test
    DL(mm): Length changes
    
    Returns:
    Depending on the operation selected,
    it returns calculated values, plots,
    advanced analysis, or saves results.
    '''
    try:
        data = pd.read_excel(file_path)
        
    except FileNotFoundError:
        print("File not found. Please check the file path.")
        return

    A0 = math.pi * (D0/2)**2

    data['stress'] = data['F (N)'] / A0
    data['strain'] = (data['DL (mm)'] - L0) / L0

    data['true_stress'] = data['F (N)'] / A0
    data['true_strain'] = np.log(1 + data['strain'])

    yield_point = data.loc[data['stress'].idxmax()]
    permanent_strain = data['strain'].iloc[-1]

    plt.figure(figsize=(12, 8))
    plt.plot(data['strain'], data['stress'], label='Engineering Stress-Strain', marker='o', color='b', linestyle='-')
    plt.plot(data['true_strain'], data['true_stress'], label='True Stress-Strain', marker='x', color='r', linestyle='--')
    plt.scatter(yield_point['strain'], yield_point['stress'], color='g', label='Yield Point')
    plt.annotate(f"Yield Point: Strain={yield_point['strain']:.2f}, Stress={yield_point['stress']:.2f}", (yield_point['strain'], yield_point['stress']), textcoords="offset points", xytext=(0,10), ha='center')
    plt.xlabel('Strain')
    plt.ylabel('Stress (MPa)')
    plt.title('Stress-Strain Curve')
    plt.grid(True)
    plt.legend()
    plt.show()

    print("Columns in the data:")
    print(data.columns)
    
    print("\nFirst few rows of the data:")
    print(data.head())

    print("\nYield Point Information:")
    print(yield_point)
    print("Permanent Strain:", permanent_strain)

file_path = 'C:/Data/Textile_Test.xlsx'  
D0 = 9
L0 = 40
Stress_Strain_Analysis(file_path, D0, L0)

import numpy as np
import matplotlib as mlp
import matplotlib.pyplot as plt
import pandas as pd


data=pd.read_excel("E:\courses\master\EVA-Gr-CB.Project\Test.Results\Mechanical.Properties\Mechanical.Properties.14\Book14.xlsx")
a=np.array(data) ##convert data to array
a1=a[11:,4]##extract elongation data
a2=a[11:,5]##extract tension data
def Stress_Strain(data,operation):
    ##drawing eva polymer tension diagram and finding max tension
    if operation=="Ultimate Tensile Stress":
        UTS=a2.max()
        return UTS
    if operation=="plot":
        font={'family':'sharif','color':'b','size':'15'}
        font2={'family':'sharif','color':'r','size':'25'}
        plt.plot(a1,a2)
        plt.title('Ethylene Vinyl Acetate',fontdict=font2)
        plt.xlabel('tension (N/mm2)',fontdict=font)
        plt.ylabel('elonfation (%)',fontdict=font)
        plt.show
  
Stress_Strain(data,'Ultimate Tensile Stress')
Stress_Strain(data,'plot')


data2=pd.read_excel("E:\courses\master\EVA-Gr-CB.Project\Test.Results\TGA\TGA1\EVA-CB\Sheet1.xlsx")
b=np.array(data2)
b1=b[1:,1]##extract temperature data
b2=b[1:,2]##extract weight data
b3=b[1:,3]##extract derivative weight data
def Thermal_Decomposition(data2,plot):
    ##drawing TGA and DTG plots to show thermal decomposition of EVA-CB composite
    if plot=="TGA":
        font={'family':'sharif','color':'b','size':'15'}
        font2={'family':'sharif','color':'r','size':'25'}
        plt.plot(b1,b2)
        plt.title('EVA-CB',fontdict=font2)
        plt.xlabel('temperature (°C)',fontdict=font)
        plt.ylabel('weight(%)',fontdict=font)
        plt.show
    if plot=="DTG":
        font={'family':'sharif','color':'b','size':'15'}
        font2={'family':'sharif','color':'r','size':'25'}
        plt.plot(b1,b3)
        plt.title('EVA-CB',fontdict=font2)
        plt.xlabel('temperature (°C)',fontdict=font)
        plt.ylabel('deriv.weight(%/°C)',fontdict=font)
        plt.show
        
Thermal_Decomposition(data2,"TGA")
Thermal_Decomposition(data2,"DTG")
        
        
        
        
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
import fsspec
import openpyxl



def Compression_Test(Data,Operator,Sample_Name,Density=0):
    """
    

    Parameters
    ----------
    Data : DataFrame
        Compression test data including stress and strain.
    Operator : str
        The action that needs to be done on the data (plot or S_max).
        plot: Plots a stress-strain diagram.
    Sample_Name : str
        Sample name.
    Density : float, optional
        Density of the sample. The default is 0.

    Returns
    -------
    float
        If the operator is given S_max, it returns maximum strength.
        If the operator is given S_max/Density, it returns specific maximum strength.

    """
    
    e=Data["e"]
    S=Data["S (Mpa)"]
    
    if Operator=="S_max":
        
        S_max=S.max()
        return S_max
        
    elif Operator=="S_max/Density" and Density!=0:
        S_max=S.max()
        S_max_Density=S_max/Density
        return S_max_Density
    
    elif Operator=="plot":
            
        font_label={   'family' : 'Times new roman' ,
                    'color' :   'black'    ,
                    'size' :   15   }

        font_legend={   'family' : 'Times new roman' ,
                     'size' :   13   }


        plt.plot(e,S,label=Sample_Name,linewidth=3)
        plt.xlabel("e",fontdict=font_label,labelpad=5)
        plt.ylabel("S (Mpa)",fontdict=font_label,labelpad=5)
    
        plt.autoscale()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.legend(frameon=False,prop=font_legend)
        plt.tick_params(axis='both', width=2)
        plt.tick_params(axis='both', which='minor', width=1)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
        plt.tick_params(axis='both', labelsize=11)
    
#testing the corroctness of the function
Data=pd.read_excel('C://Users//Asus//Desktop//University//Master of Material Science//Master Project//Analysis Results//Compression//Test 26//Excel Files//G20525-4-1.xlsx')
Data1=pd.read_excel('C://Users//Asus//Desktop//University//Master of Material Science//Master Project//Analysis Results//Compression//Test 26//Excel Files//PBS20525-1-1.xlsx')
Compression_Test(Data,"plot","G20525-4")
Compression_Test(Data1,"plot","PBS20525-1")
Compression_Test(Data,"S_max","G20525-4")
Compression_Test(Data,"S_max/Density","G20525-4",0.61)






def DMTA_Test(Data,Operator,Sample_Name):
    """
    

    Parameters
    ----------
    Data : DataFrame
        DMTA test data including storage modulus, loss modulus and tanδ.
    Operator : str
        The action that needs to be done on the data (storage_max, loss_max, tan_max, plot_storage, plot_loss or plot_tan).
    Sample_Name : str
        Sample name.

    Returns
    -------
    float
        If the operator is given storage_max, it returns maximum storage modulus.
        If the operator is given loss_max, it returns maximum loss modulus.
        If the operator is given Tan_max, it returns maximum Tanδ.

    """
    
    Frequency=Data2["Frequency (Hz)"]
    Storage_Modulus=Data2["E'-Storage Modulus (Mpa)"]
    Loss_Modulus=Data2.iloc[:,13].copy()
    Tanδ=Data2["Tanδ"]

    if Operator=="storage_max":
        Storage_max=Storage_Modulus.max()
        return Storage_max
    
    elif Operator=="loss_max":
        Loss_max=Loss_Modulus.max()
        return Loss_max
    
    elif Operator=="tan_max":
        Tan_max=Tanδ.max()
        return Tan_max
    
    elif Operator=="plot_storage":
        
        font_label={   'family' : 'Times new roman' ,
                    'color' :   'black'    ,
                    'size' :   15   }

        font_legend={   'family' : 'Times new roman' ,
                     'size' :   13   }


        plt.plot(Frequency,Storage_Modulus,label=Sample_Name,linewidth=3)
        plt.xlabel("Frequency (Hz)",fontdict=font_label,labelpad=5)
        plt.ylabel("E'-Storage Modulus (Mpa)",fontdict=font_label,labelpad=5)
    
        plt.autoscale()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.legend(frameon=False,prop=font_legend)
        plt.tick_params(axis='both', width=2)
        plt.tick_params(axis='both', which='minor', width=1)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
        plt.tick_params(axis='both', labelsize=11)
        
    elif Operator=="plot_loss":
        
        font_label={   'family' : 'Times new roman' ,
                    'color' :   'black'    ,
                    'size' :   15   }

        font_legend={   'family' : 'Times new roman' ,
                     'size' :   13   }


        plt.plot(Frequency,Loss_Modulus,label=Sample_Name,linewidth=3)
        plt.xlabel("Frequency (Hz)",fontdict=font_label,labelpad=5)
        plt.ylabel("E''-Loss Modulus (Mpa)",fontdict=font_label,labelpad=5)
    
        plt.autoscale()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.legend(frameon=False,prop=font_legend)
        plt.tick_params(axis='both', width=2)
        plt.tick_params(axis='both', which='minor', width=1)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
        plt.tick_params(axis='both', labelsize=11)    
    
    elif Operator=="plot_tan":
        
        font_label={   'family' : 'Times new roman' ,
                    'color' :   'black'    ,
                    'size' :   15   }

        font_legend={   'family' : 'Times new roman' ,
                     'size' :   13   }


        plt.plot(Frequency,Tanδ,label=Sample_Name,linewidth=3)
        plt.xlabel("Frequency (Hz)",fontdict=font_label,labelpad=5)
        plt.ylabel("Tanδ",fontdict=font_label,labelpad=5)
    
        plt.autoscale()
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
        plt.legend(frameon=False,prop=font_legend)
        plt.tick_params(axis='both', width=2)
        plt.tick_params(axis='both', which='minor', width=1)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
        plt.tick_params(axis='both', labelsize=11)


#testing the corroctness of the function
Data2=pd.read_excel('C://Users//Asus//Desktop//University//Master of Material Science//Master Project//Analysis Results//DMTA//Compression Mode//Test 1//Excel Files//G20525-4.xlsx')
DMTA_Test(Data2,"storage_max","G20525-4")
DMTA_Test(Data2,"loss_max","G20525-4")
DMTA_Test(Data2,"tan_max","G20525-4")
DMTA_Test(Data2,"plot_storage","G20525-4")
DMTA_Test(Data2,"plot_loss","G20525-4")
DMTA_Test(Data2,"plot_tan","G20525-4")


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#Please input the data file in the line below!
data=pd.read_excel('C:/Users/Lenovo/Downloads/ai_2024.xlsx')



#First_Function

def Find_Max_Vertical_Velocity(data):
    '''
    This function gets the data results of flow velocity simulation for natural convection of a molten in a cavity.
    and as the output, it shows the plot of flow velocity vs x
    and also it returns the maximum velocity and the location where that velocity belogs to.
    '''
    
    x=np.array(data['x(m)'])
    u=np.array(data['u(m/s)'])
    u_max=np.max(u)#the max value of the velocity
    index_max=np.argmax(u)#the index of which this value exists in
    loc_max=x[index_max]#the location of the maximum velocity
    print('The maximum value of Flow Velocity for this problem is:',u_max)
    print('Also this maximum value occurs in this location:',loc_max)
   
    plt.scatter(x,u)
    plt.title('Flow Velocity',c='blue',family='times new roman',size=20,pad=20)
    plt.xlabel('x (m)',size=15,c='green')
    plt.ylabel('u velocity (m/s)',size=15,c='green')
    plt.show()
    
    return u_max,loc_max


Find_Max_Vertical_Velocity(data)
    
    
#Second_Function
        
def Solidification_Start(data,Temp_sol):
    '''
    This function gets the temperature values of a molten along with its solidus temperature as inputs,
    and if the solidification process has started it returns true and if not it returns false.
    it also plots the temperature values in the center line.
    note that this function has simplified the problem and many other parameters should be taken into consideration.
    '''
    
    x=np.array(data['x(m)'])
    y=np.array(data['T(K)'])
    plt.scatter(x,y)
    plt.title('Temperature Profile',c='blue',family='times new roman',size=20,pad=20)
    plt.xlabel('x (m)',size=15,c='green')
    plt.ylabel('Temperature (K)',size=15,c='green')
    plt.show()
    Temp_Min=np.min(y)
    if Temp_Min>Temp_sol:
        print('The solidification process has not started yet.')
    else:
        print('The solidification process has started.') 
        
    
    if Temp_Min>Temp_sol:
        return False
    else:
        return True
    

Solidification_Start(data, 466)  
    

import numpy as np
import pandas as pd
import openpyxl
import matplotlib as mlp
import matplotlib.pyplot as plt
    

#===========
font_a ={ 'family': 'serif',
             'color' : 'red',
             'size' : 14
        }  

#===========


'''
===========================================First

تابع توزیع احتمال دو جمله ای
برای حساب کرن توزیع احتخمال یک پیشامد
این تابع برای محاسبه احتمال پیشامد بصورت یک و صفر است
مثلا برای حساب کردن احتمال اینکه با 10 با پرتاب سکه چند بار 6 ببینیم
که یک صورت شکست دارد و یک صورت پیروزی مثلا شیر پیروزی و خط شکست باشد
_______________________

فرمول محاسبه تابه احتمال دو جمله ای
P(x):(tarkib X az N)*(p**x)*(q**n-x)
tarkib X az N= n! / x! (n-x)!

------> P(x): احتمال وقوع پیشامد
------> n = تعداد تکرار آزمایش
------> p = 1 piroozi
------> q = p-1 faild
------> x = moteghayer tasadofi

توابع مورد نیاز 
Factoriel
'''

def fact(x):
    if x==0:
        return 1
    else:
        return x * fact(x-1)
'''    
>>>>>>>>>>>>tabe avalie
def Bino_m(x,n,p=0.5):
    q = 1 - p
    fact_n = float(fact(n))
    fact_x = float(fact(x))
    fact_nx= float(fact(n-x))
    tarkib_x_n = fact_n / (fact_x * fact_nx)
    px = tarkib_x_n * (p**x) * ((q)**(n-x))
    
    
    #plt.plot(args, kwargs)
    plt.show()
    return px
'''




def Bino_m_rasm(x,n,p=0.5):
    '''
        Parameters
    ----------
    x : int
        addad mored azmayesh.
    n : int
        tedad tekrar azmayesh
    p : float optional
      ehtemal pirouzi. The default is 0.5.

    Returns
    -------
    plot
        neshan dadan bordar va natije
    float
        mizan ehtemal pirouzi.

    '''

    first_x=x
    q = 1 - p
    x=np.arange(1,n+1)
    tarkib_x_n = np.array([])
    fact_n = float(fact(n))
    px=np.array([])
    a=np.array(())
    for i in x:
        fact_x = float(fact(i))
        n2=n-i
        fact_nx= float(fact(n2))
        tarkib_x_n = float (fact_n / (fact_x * fact_nx))
        px_x = float( tarkib_x_n * (p**i) * (q**n2) )
        px=np.append(px,px_x)
        a=np.append(a,i)
        if i==first_x:
            plt.plot(a[first_x-1], px[first_x-1],marker="o",mec='g')
        else:
            plt.plot(a,px)
                
    return px[first_x-1],plt.show()

#================test:
v=Bino_m_rasm(6,10)

       
"""
===========================================second
mizan sakhti ab
عوامل موثر بر سختی آب میزان
در درجه اول میزان کلسیوم کربنات و منیزیوم در آن است
و فرمول آن به شرح زیر است
ppm(sakhti ab-->mg/L) = (Mg * 4.12) + (Ca * 2.5)
در صورتی که عدد به دست آمده از 60 کمتر باشد آب نرم
اگر بین 60 تا 120 باشد آب متوسط
بین 120 تا 180 باشد آب سخت
و بیشتر از 180 آب بسیار سخت میباشد
همچنین میازن مس آب باید کمتر از 20 ملیگ گرم در لیتر - نیکل کمتر از 10 زینک کمتر از 10 
آرتوفسفات کمتر از 100 و سیانید کمتر از 2 باشد
در غیر اینصورت آب قابل استفاده نیست.
====> Mg manyazio,
====> ca calsium carbonat 
در این تابع ورودی از اکسل دریافت و نتیجه آزمایش ها رسم میشود 
در صورتیکه آب غابل استفاده نباشد در لیست مجزایی دی خروجی نمایشگ داده خوادهد شد
"""

data = pd.read_excel('C:\\Users\\PARVANE-PC\\Desktop\\samples.xlsx')



def Sakhti_ab( data):
    '''
    Parameters
    ----------
    data : exel
        شامل تمامی داده های مورد نیاز برای بررسی سختی آب میباشد.

    Returns
    -------
    data : plot
        plot mizan sakhti ab dar nemoone ha.
    arraye
        jadval nahayi.

    '''
    
    unavalable_water=[]

    drop_index= data[data['Cu']>20].index
    for i in drop_index:
        unavalable_water.append(data.loc[[i],['name']])
    data= data.drop(drop_index)
    
    drop_index= data[data['Ni']>10].index
    for i in drop_index:
        unavalable_water.append(data.loc[[i],['name']])
    data= data.drop(drop_index)
    
    drop_index= data[data['Zn']>10].index
    for i in drop_index:
        unavalable_water.append(data.loc[[i],['name']])
    data= data.drop(drop_index)
    
    drop_index= data[data['pyro']>100].index
    for i in drop_index:
        unavalable_water.append(data.loc[[i],['name']])
    data= data.drop(drop_index)
    
    drop_index= data[data['Cya']>2].index
    for i in drop_index:
        unavalable_water.append(data.loc[[i],['name']])
    data= data.drop(drop_index)
    #======
    
    mg = np.array(data['Mg'])
    ca = np.array(data['Ca'])
    names= np.array(data['name'])

    ppm = np.array([])
    ca = (ca*2.5)
    mg = (mg*4.12)
    ppm = ca + mg
    plt.bar(ppm,names)
    
    return data,plt.show()

#================= test

new1 =Sakhti_ab(data)

'''
در این تابع با دسترسی توابع پانداس تونستم به متغیرهای داخل دیتا فرم دسترسی پیدا کنم
مقادیر که کلا در ازمایش اولیه حذف میشه رو از دیتا فرم حذف کردم
و درجه سختی رو با جدول باقی مانده انجام دادم
میخاستم جدول نهایی خروجی ای از میزان سختی باشه وابسته به شهر ولی نمیتونستم مقدار آرایه ppm  به حالت بازه تغییر بدم 
یعنی نتونستم دسترسی به اندیس پیدا کنم
و با حلقه for دسته بندی کنم


'''

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd


# 1-----------------------------------------------------------------------------
w=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight CrN S.xlsx')


'''
bayad adresse mahalle zakhireye dadeha dar w vared shavad.
'''


def Wear_Rate(w,work,S,F):
    '''
    w bayad shamele 2 soton bashad ke dar yek soton vazn nemoone ghabl az test andazegiri shavad ,
    dar yek soton vazne nemone pas az etmame test.
    
    S haman masate laghzesh hine azmone sayesh ast. S bayad bar hasbe metr bashad.
    
    F barabare niroyee ast ke be pin vared shode ast , azmon ba an anjam shode ast. F bayad bar hasbe newton bashad
    '''
    w.columns
    wb=np.array(w['weight befor test'])
    # khate bala soton avval ra joda mikonad yani vazne nemone ghabl az azmoon.
    wa=np.array(w['weight after test'])
    # khate bala soton dovvom ra joda mikonad yani vazne nemone baa'd az azmoon.
    wl=np.subtract(wb,wa)
    # khate bala kaheshe vazne nemoone pas az azmoon ra hesab mikonad.
    m=wl.mean()
    # khate bala miangine kaheshe vazne nemoone ra hesab mikonad.
    if work=='wear rate':
        WR= m/(S*F)
        # WR yani Wear Rate va nerkhe sayesh ra hesab mikonad.
        return WR
    
    
a=Wear_Rate(w,'wear rate',300,5)
print('nerkhe sayeshe shoma barabar ast ba:', a)




w=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=0)


def Wear_Rate(w,work,S,F):
        '''
        w bayad shamele 2 soton bashad ke dar yek soton vazn nemoone ghabl az test andazegiri shavad ,
        dar yek soton vazne nemone pas az etmame test.
        S haman masate laghzesh hine azmone sayesh ast. S bayad bar hasbe metr bashad.
        F barabare niroyee ast ke be pin vared shode ast , azmon ba an anjam shode ast. F bayad bar hasbe newton bashad
        '''
        w.columns
        wb=np.array(w['weight befor test'])
        # khate bala soton avval ra joda mikonad yani vazne nemone ghabl az azmoon.
        wa=np.array(w['weight after test'])
        # khate bala soton dovvom ra joda mikonad yani vazne nemone baa'd az azmoon.
        wl=np.subtract(wb,wa)
        # khate bala kaheshe vazne nemoone pas az azmoon ra hesab mikonad.
        m=wl.mean()
        # khate bala miangine kaheshe vazne nemoone ra hesab mikonad.
        if work=='wear rate':
            WR= m/(S*F)
            # WR yani Wear Rate va nerkhe sayesh ra hesab mikonad.
            return WR






w0=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=0)
w1=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=1)
w2=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=2)
w3=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=3)
w4=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=4) 
#sheet name har safhe az excel ra baz mikonad.

list_1=[w0,w1,w2,w3,w4]
# yek list az dadeha misazim.
def Wear_Bar(list_1,work='bar'):
    # manzoor az bar nemoodare mile ee ast.
    new_list=[]
    if work=='bar':
        for i in list_1:
            # be komake halgheye for nerkhe sayeshe tamame nemoone ha ra hesab mikonim.
            aa=Wear_Rate(i,'wear rate',300,5)
            new_list.append(aa)
            # nerkhe sayeshe har nemoon ra be yek lidt ezafe mikonim.
    x=np.array(['A','B','C','D','E'])
    y=np.array(new_list)
    plt.bar(x,y,color='g',width=0.5)
    # nemidare mile ee ra ba dastoore plt.bar rasm mikonim.
    plt.title('The amount of wear rate for samples',c='g',size=14)
    plt.ylabel('Wear Rate(mg/N.m',size=12,color='k')
    plt.show()
    
    
w0=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=0)
w1=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=1)
w2=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=2)
w3=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=3)
w4=pd.read_excel('/Users/intel/Desktop/poject 2/wear weight.xlsx',sheet_name=4) 

list_1=[w0,w1,w2,w3,w4]
Wear_Bar(list_1,'bar')









# 2-----------------------------------------------------------------------------
font_1={'family': 'serif',
        'color': 'b',
        'size':14}

font_2={'family': 'serif',
        'color': 'k',
        'size':12}
# 2 ta font misazim ke betavanim estefade konim.

b=pd.read_excel('/Users/intel/Desktop/poject 2/polarization substrate.xlsx')
def Polarization(b,work):
    '''
    d: dadehaye azmayeshgah shmele ghegaliye jaryan va potansiel hastand.
    '''
    b.columns
    cd=np.array(b['Current density'])
    # cd haman ghegaliye jaryan ast.
    cp=np.array(b['Potential'])
    # cp haman potansiel ast.
    cdl=np.log(cd)
    # dadehaye ghegaliye jaryan bayad logaritni shavand.
    if work=='plot':
        # in ghesmat marboot be rasme nemoodare polarizasion ast.
        plt.plot(cdl,cp,c='r',linewidth=3)
        plt.title('Polarization Curve',fontdict=font_1,pad=15)
        plt.xlabel('Log i(A/cm2)',fontdict=font_2)
        plt.ylabel('E(V Ag/AgCl)',fontdict=font_2)
        # dar 3 khate bala az font haye sakhte shode estefade shod.
        plt.show()
    if work=='corrosion potential':
        # in ghesmat marboot be mohasebeye potansiele khordegi ast.
        c=cdl.min()
        # meghdare minimume ghegalie jaryan ra miyabim ta potansiele mirboot be an ra peyda konim ke an,
        # potansiele khordegi ast.
        d=np.where(cdl==c)
        # baraye yaftane mahale minimume ghegalie jaryan.
        e=b.loc[d]
        f=np.array(e)
        g=f[0,1]
        return g

h=Polarization(b, 'corrosion potential')
print('potansiele khordegie shoma barabar ast ba:',h)


Polarization(pd.read_excel('/Users/intel/Desktop/poject 2/polarization substrate.xlsx'), 'plot')
# ya
# Polarization(b, 'plot')
        
        



         


import numpy as np#for array and calculation
import pandas as pd#for data
import matplotlib.pyplot as plt#for plot

data_path = r'C:\Users\DataSystem.ir\Desktop\Stress_strain.xlsx'#here we get the path and file name
data = pd.read_excel(data_path)#read the data

#in our excel file we have Forza standard(N) and deformation and area is constantant which is 5/29916E-05
#first we write def get force and area to calculate stress

area = 5 / 29916E-05 # Constant area
#constant number for our calculation so I define it here

def Stress_Cal(data, area):
    N_data = data['N']#access data in column N which is related to force
    x = np.array(N_data / area / 1000000)#creat array of them
    # the formula for stress= Force/area and it has division by 1000000 for Mpa
    return x

def Strain_Cal(data):#y is our strain and it divided by 100 since it is percent of deformation
    Deformation_data = data['Deformation']
    y = np.array(Deformation_data / 100)
    return y

def Stress_Strain(x, y, operation):#calculate modulus
    m = np.divide(x, y, out=np.zeros_like(x), where=y != 0) # in this part first I wrote if y==o continue but i search in the web and they recommended this id strain is zero continue and just remain it zerp
    print(m)#calculate modulus by using divivision in array

#in this part which is related to operation we calculate min, max and plot the diagram

    if operation == 'min':#calulate min stress and strain
        min_stress = x.min()
        min_strain = y.min()
        print('Min strain=',min_strain,'Min stress=',min_stress)
        
    elif operation == 'max':#calculate max strain and straee
        max_stress = x.max()
        max_strain = y.max()
        print('Max strain=',max_strain,'Max stress=', max_stress)
        
    elif operation == 'plot':# draw plot
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, 'o-', label='Stress-Strain Relationship',c='g')
        plt.xlabel('Strain')
        plt.ylabel('Stress')
        plt.title('Stress-Strain Diagram')
        plt.legend()#for better understangin when the data are alot and close to each other
        plt.grid(True)
        plt.show()
        
        
def Slope_Cal(x, y):# in this part I search in the web
    # Filter out the zero strain values to avoid division by zero because if it divided by zero in would be INF
    fil_indices = y != 0 # just why becase just strain is impo for us and if strain would be zero can cause the problem for us
    fil_x = x[fil_indices]
    fil_y = y[fil_indices]

    # Fit a linear model (polynomial of degree 1) and it should be one because the slope in linear part is impo
    slope, intercept = np.polyfit(fil_y, fil_x, 1)
    return slope, intercept

operation = input("Enter max, min, or plot: ")   #I suppose that user enter the operator   
y = Strain_Cal(data)
x = Stress_Cal(data, area) 
# Then, in your main code block, call this function with your data:  
slope, intercept = Slope_Cal(x, y)
print("Slope (Modulus of Elasticity):", slope)  
Stress_Strain(x, y, operation)




#the slope in my data are not correect because I limit my data and I just use first 100 data and the linear part for calculate slope may be in the the after these

#------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy import stats
#------------------------------------------------------------------------------
# import data
xrd = pd.read_excel(r'D:\6. ai in mse\library\quiz\xrd.xlsx')
# data description
xrd .columns
#------------------------------------------------------------------------------
# font
font_title = { 'family' : 'serif','color' :'black' , 'size' : 15 }
font_x = {'family' :  'serif','color' :   'black', 'size' : 10 }
font_y = { 'family' : 'serif','color' :   'black', 'size' : 10}    
#------------------------------------------------------------------------------
#def function
def Xrd_Analysis(data, operation):
    
    '''
    XRD is a nondestructive technique that provides detailed
    information about the crystallographic structure, chemical 
    -----------------------------------------------------------
    composition, and physical properties of a material
    this function return max intensity of xrd pattern and degree 
    of difraction. also based on user requierment return specified plot.
    statistics analysis is presented too. 
    ----------------------------------------------------------
    Parameters
    ----------
    data : DataFrame
        data set include tow features:
            1. degree (2* degree of difraction)
            2. intensity
    operation : str
        determine an operation which apply on data
        possible values:
            information
            max intensity
            scatter plot
            line graph
    Returns
    -------
    max intensity and degree of max intensity 
    or
    plot

    '''


    if operation == 'max intensity':
        max_intensity = data['intensity'].max()
        specified_intensity =  data['intensity'] == max_intensity
        related_vlues = data [specified_intensity]
        return related_vlues
    
    elif operation == 'scatter plot':
        y = np.array(data['intensity'])
        x = np.array(data['degree'])
        plt.scatter(x, y, marker = '^', color = 'red', edgecolors='brown')
        plt.title('XRD pattern',fontdict=font_title )
        plt.xlabel('Degree',fontdict=font_x)
        plt.ylabel('Intensity',fontdict=font_y)
        plt.show()
    
    elif operation == 'line graph':
        y = np.array(data['intensity'])
        x = np.array(data['degree'])
        plt.plot(x,y, color = 'green', linewidth = 0.5)
        plt.title('XRD pattern',fontdict=font_title )
        plt.xlabel('Degree',fontdict=font_x)
        plt.ylabel('Intensity',fontdict=font_y)
        plt.show()
        
    elif operation == 'line graph fill between':
        y = np.array(data['intensity'])
        x = np.array(data['degree'])
        plt.fill_between( x, y, color="#C8D700" , alpha = 0.3) 
        plt.plot(x, y, color='#36BD00', linewidth = 0.5)
        plt.title('XRD pattern',fontdict=font_title )
        plt.xlabel('Degree',fontdict=font_x)
        plt.ylabel('Intensity',fontdict=font_y)
        plt.show()

#------------------------------------------------------------------------------
Xrd_Analysis(xrd, 'max intensity')   
Xrd_Analysis(xrd, 'scatter plot')   
Xrd_Analysis(xrd, 'line graph')   
Xrd_Analysis(xrd, 'line graph fill between')   
#------------------------------------------------------------------------------
#import data
cof = pd.read_excel(r'D:\6. ai in mse\library\quiz\cof.xlsx')
cof.columns
#------------------------------------------------------------------------------
# def function
def Statitical_Analysis(data, operation):
    '''
    this function calculate quantile, IQR, min, max, median,and
    zscore for each features of data set. also it is presented
    plot.
    -----------------------------------------------------------
    Parameters
    ----------
    data : data frame
        
    opertaion : str
        possible values:.
        1. statistics
        2. histogram
        3. correlation
        4.pairplot
    Returns
    -------
    1. quantile
    2. min
    3. max
    4. median
    5. zscore
    6. determined plot
    '''
    
    if operation == 'statistics':
       for c in data.columns:
           if (data[c].dtype == int) | (data[c].dtypes == float):

                Q1 = data[c].quantile(.25)
                Q3 = data[c].quantile(.75)
                IQR = Q3 -Q1
                
                min_value = data[c].min()
                max_value = data[c].max()
                median_value = data[c].median()
                
                z_score = np.abs(stats.zscore(data[c]))
                
                print('feature name : ', c,
                      '\n min = ', min_value, 
                     'max = ',max_value ,
                     'median = ',median_value ,
                     'Q1 = ', Q1, 'Q3 = ', Q3, 'IQR = ', IQR,
                     'ZScore = ', z_score)
    
    for c in data.columns:
        if operation == 'histogram':
            
            if (data[c].dtype == int) | (data[c].dtypes == float):
                
                plt.hist(data[c], label = c, color = 'green',
                         edgecolor = 'black', linewidth = 0.5)

                plt. legend()
                plt.title('distribution', fontdict=font_title)
                plt.xlabel('bins', fontdict = font_x)
                plt.ylabel('frequency', fontdict = font_y)
                plt.show()
    
    if operation == 'correaltion':
        plt.figure(figsize=(18,12), dpi= 200)

        sns.heatmap(data.corr(), xticklabels=data.columns, 
            yticklabels=data.columns, center=0, annot=True,cmap = 'coolwarm')
                
        plt.title('correalation', fontdict=font_title)
        plt.show()
    if operation == 'pairplot':
         plt.figure(figsize=(18,12), dpi= 200)
         sns.pairplot(data, vars = data.columns, markers = '*', kind = 'reg')
         plt.title('relation between columns', fontdict=font_title)
         plt.show()
         

#------------------------------------------------------------------------------
Statitical_Analysis(cof, 'statistics')     
Statitical_Analysis(cof, 'histogram')    
Statitical_Analysis(cof, 'correaltion')  
Statitical_Analysis(cof, 'pairplot')  
#------------------------------------------------------------------------------


import pandas as pd
import matplotlib.pyplot as plt
data1=pd.read_excel('C:/Users/Gandi/Desktop/Boo.xlsx')
data2=pd.read_excel('C:/Users/Gandi/Desktop/Bo.xlsx')

def Blood_Pressure(data1,operation):
    '''
 Blood is divided into two main values: systolic blood pressure and diastolic blood pressure.
 Systolic blood pressure represents the blood pressure at the moment the heart muscle contracts,
 while diastolic blood pressure represents the blood pressure at the moment the heart muscle relaxes.

 This function gives us their average systolic and diastolic blood pressure
 And it shows us the blood pressure chart of 40-year-old Balinese people
    Parameters
    ----------
    data1 : int
        systolic,diastolic
    operation : strig
    operator
    Returns
    -------
    None.

    '''
    a=data1['Systolic']
    b=data1['Diastolic']
    
   
    if operation=='average1':
     c=b.mean()
     return c
    if operation=='average2':
        c=a.mean()
        return c
    if operation=='plot':
        plt.plot(a,b,marker='o',ms=10,mec='r')
        plt.xlabel('Systolic')
        plt.ylabel('Diastolic')
        plt.title('Blood_Pressure')
        plt.show()
        
def Pulse(data2,operation):
    '''
    This function gives us the maximum and minimum pulse of people over forty years old 
    in the test and draws the pulse graph of people.
    And by using the frequency of the pulses, it is obtained

    Parameters
    ----------
    data2 : int
        Pulse
    operation : string
        operator

    Returns
    -------
    None.

    '''
    p=data2['pulse'] 
    if operation=='maximum':
        c=p.max()
        return c
    if operation=='minimum':
        c=p.min()
        return c
    if operation=='plot':
        plt.hist(p)
        plt.title('Pulse')
        plt.show()
        

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits import mplot3d
import math
import statistics

filee=pd.read_excel('C:/Users/HP/Desktop/res.xlsx')
filee1= pd.read_excel('C:/Users/HP/Desktop/dls.xlsx')

#first function:

def Color_Feature(filee,kar):
    
    '''

    Parameters
    ----------
    filee : DataFrame
        data from excel.
    kar : str
        plot or calculate?

    Returns
    a number as a delta_E : array or plot a diagram

    '''
    l=filee['l']
    a= filee['a']
    b= filee['b']
    if kar=='calculate':
        l_=[]
        a_=[]
        b_=[]
        for i in l:
            l_.append((i-l[5])**2)
        for i in a:
            a_.append((i-a[5])**2)
        for i in b:
            b_.append((i-b[5])**2)
        
        sum1 = np.add(l_,a_)
        sum_=np.add(sum1,b_)
        squrt_=[]
        for i in sum_:
            squrt_.append(math.sqrt(i))
        delta_E= np.array(squrt_)
        return delta_E
    if kar=='plot':
         fig = plt.figure()
         ax= fig.add_subplot(111,projection='3d')
         ax.set_xlabel('l*')
         ax.set_ylabel('a*')
         img= ax.scatter(l,a,b,color='r')
         ax.set_title('l* a* b* of color')
Color_Feature(filee,'plot')

#second function:

def Particles_Size(filee1,kar1):
    
    '''
    

    Parameters
    ----------
    filee1 : DataFrame
        Data from excel.
    kar1 : str
        plot or calculation.

    Returns
    a number as an average size or plot a diagram.
    

    '''
    x= filee1['size']
    y= filee1['distribution']
    
    if kar1=='calculate':
        size=filee1['size']
        average=statistics.mean(size)
        print('your average size is:')
        return average
    if kar1=='plot':
        x= filee1['size']
        y= filee1['distribution']
        font_x={'family':'Times New Roman',
             'color':'k',
              'size':'15'}
        font_y={'family':'Times New Roman',
             'color':'k',
              'size':'15'}
        plt.xlabel('size(nm)',fontdict=font_x)
        plt.ylabel('intensity(%)',fontdict=font_y)
        plt.plot(x,y)
        
    
Particles_Size(filee1,'plot')  

    

def Price_Change(data,operation):
    '''Parameters
    ----------
    data : Csv File
        data ha darbare vizhegi haye gheymat yek saham dar yek rooz mibashad
    descrription
     
    operation : az beyn "mohasebe" va "plot"yeki ra entekhab konid
    in function bishtarin va kamtarin gheymat saham dar yek rooz ra migirad
    va damane taghyirat gheymat ra ba "plot"rasm mikonad va ba "mohasebe"
    return mikonad
    '''
    if operation=="plot":
         plt.plot(data.High - data.Low)
    if operation=="mohasebe":
        data["taghyir"]=data['High']-data['Low']
        new_data=data[["High","Low","taghyir"]]
        return new_data["taghyir"] 
       
df=pd.read_csv("E:/codebasis/pandas/july2017apple.csv")         
Price_Change(df,"mohasebe")
Price_Change(df,"plot")
def New_Case_Corona_Propotion(data,operation):
    ''' Parameters
    ----------
    data : Csv File
        data hayi darbare amar corona virus dar keshvar haye mokhtalef
    operation : az beyn "plot" va "nesbat" yeki ra mitavanid entekhab konid
        in tabe nesbat new case haye coronavirus be case haye ghabli ra 
        ba "plot" rasm mikonad va ba "nesbat" return mikonad

    '''
    if operation=="plot":
        plt.plot(data.Cases-data.New_cases)
    if operation=="nesbat":
        data["New_Propotion"]=data["New_cases"]/data["Cases"]
        return data.New_Propotion
df=pd.read_csv("C:/Users/MashadService.ir/Desktop/corona_march2020_mlcourse.csv",index_col="Country") 
New_Case_Corona_Propotion(df,"plot") 
New_Case_Corona_Propotion(df,"nesbat")     
        




import numpy as np
import pandas as pd
import matplotlib as mlp
import matplotlib.pyplot as plt
def Load_Position_Convertor(Data,Operation,Area,Length):
    '''
    This function receives an input file containing Load-Position data, as well as the cross-sectional area and gauge length of the part, and according to the user's needs, it can:
    1-Load-Position Curve (LPC)
    2-Stress-Strain Calculation (SSCal)
    3-Stress-Strain Curve (SSC)
    4-Normalized Stress-Strain Calculation (NSSCal)  
    5-Normalized Stress-Strain Curve (NSSC) 
    6-Energy Absorption Density Calculation (EADCal)

    Parameters
    ----------
    Data : xlsx
        Need two columns containing Load-Position information.
    Operation : str
        It specifies the process that should be performed on the entered information.
    Area : float
        The surface examined in the tensile or compression test.
    Length : TYPE
        Gauge length checked in tension or compression test.

    Returns
    -------
    EAD : float
        Energy absorption desity of meta-materials such as metal foams.

    '''

    Data=pd.read_excel('C://Users//Dr Computer//Desktop//Data For Project 2.xlsx')
    Load=np.array(Data['Load (kN)'])
    Position=np.array(Data['Position (mm)'])

    
    if Operation=='Load-Position Curve' or 'LPC':
        plt.plot(Position,Load,c='teal',lw='1')
        title_font={'family':'times new roman','color':'black','size':'14'}
        label_font={'family':'times new roman','color':'black','size':'12'}
        plt.title('Load-Position',fontdict=title_font,loc='center',pad=10)
        plt.xlabel('Position (mm)',fontdict=label_font,labelpad=5)
        plt.xlim(0,np.max(Position))
        plt.ylabel('Load (kN)',fontdict=label_font,labelpad=5)
        plt.ylim(0,np.max(Load))
        plt.grid(linewidth=0.5,color='grey',alpha=0.4)
        plt.show()
    elif Operation=='Stress-Strain Calculation' or 'SSCal':
        Strain=Position/Length
        Stress=(Load*1000)/Area
        Stress_Strain=np.dstack((Strain,Stress))
        return Stress_Strain
    elif Operation=='Stress-Strain Curve' or 'SSC':
        Strain=Position/Length
        Stress=(Load*1000)/Area
        plt.plot(Strain,Stress,c='teal',lw='1')
        plt.title('Stress-Strain',fontdict=title_font,loc='center',pad=10)
        plt.xlabel('Strain (-)',fontdict=label_font,labelpad=5)
        plt.ylabel('Stress (MPa)',fontdict=label_font,labelpad=5)
        plt.grid(linewidth=0.5,color='grey',alpha=0.4)
        plt.show()
    elif Operation=='Normal Stress-Strain Calculation' or 'NSSCal':
        N_Strain=Strain/Strain.max()
        N_Stress=Stress/Stress.max()
        N_Stress_Strain=np.dstack(N_Strain,N_Stress)
        return N_Stress_Strain
    elif Operation=='Normal Stress-Strain Curve' or "NSSC":
        N_Strain=Strain/Strain.max()
        N_Stress=Stress/Stress.max()
        plt.plot(N_Strain,N_Stress,c='teal',lw='1')
        plt.title('Normal Stress-Strain',fontdict=title_font,loc='center',pad=10)
        plt.xlabel('Normal Strain (-)',fontdict=label_font,labelpad=5)
        plt.xlim(0,1)
        plt.ylabel('Normal Stress (-)',fontdict=label_font,labelpad=5)
        plt.ylim(0,1)
        plt.grid(linewidth=0.5,color='grey',alpha=0.4)
        plt.show()
    elif Operation=='EAD Calcultion' or 'EADCal':
        EAD=np.trapz(Stress,Strain)
        return EAD


        
Load_Position_Convertor( '','EAD Calculation', 10, 10) 

#1------------------------------------------------


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mlp
import math
 



 
def Nemodar_Sahmi(x,y):
    x=np.arange(0, 2*(np.pi), 0.1) 
    y=np.sin(x)
  
        
plt.plot(x,y,marker='o',ms=5,mec='tan',mfc='black',ls='--',c='tan',lw='2') 

font_title={ 'color':'olive',
             'size':25}

font_xlabel={'color':'black',
              'size':15}

font_ylabel={'color':'black',
              'size':15}

plt.title('Nemodar Sahmi',font_title)

plt.xlabel( 'X',font_xlabel)
plt.ylabel( 'Y', font_ylabel)
plt.grid()
plt.show()
  

 


#2--------------------------------------------------



def Quadratic_Function(a,b,c) :
    x=np.arange(-1000,1000,10 )
    a=2
    b=15
    c=10
    y=a*(x**2)+b*x+c
       
plt.plot(x,y,marker='x',ms=10, ls='-.',color='r',mfc='r',mec='g',lw='1')
plt.title("Quadratic Function",color='r',pad=20)
plt.xlabel("Values of x",color='b')
plt.ylabel("Values of y",color='b')
plt.show()






#3--------------------------------------------------------




plt.plot(np.sin(x),label='sin(x)',color='b')
plt.plot(np.cos(x),label='cos(x)',color='r')
plt.legend()
plt.title('Sin_Cos',pad=20,loc='left')
plt.grid()
plt.show()

    
  


 















f_loc='//c://Users//Rayan//Desktop//'
data=pd.read_excel(f_loc)










def Brain_Data():
    """
    This function prompts the user to input values for sudiom and potasiom along with time increments,
    creates a DataFrame with the provided data, and returns it.
    """
    # Initialize an empty dictionary to store the data
    data = {
          'sudiom':[],
          'potasiom':[],
          'time': []  
    }

    # Prompt the user to enter the number of entries
    Num_Entries = int(input('Your entries ion ? --> '))
    time = 0  # Initial time value
    counter = 1  # Initial counter value

    # Loop to input data from the user
    for i in range(Num_Entries):
        sudiom = float(input(f'sudiom {counter}: '))  # Prompt for sudiom input
        potasiom = float(input(f'potasiom {counter}: '))  # Prompt for potasiom input
        counter += 1  # Increment counter

        # Append the input values to the data dictionary
        data['sudiom'].append(sudiom)
        data['potasiom'].append(potasiom)
        data['time'].append(time)  # Append the time value
        time += 0.1  # Increment time by 0.1 for each entry

    # Convert the data dictionary to a pandas DataFrame and return it
    df = pd.DataFrame(data)
    return df

def Brain_Plot_3D(df):
    """
    This function takes a DataFrame containing sudiom, potasiom, and time data and plots
    them in both 2D and 3D. It also compares the mean values of sudiom and potasiom to decide
    which one to plot in 3D.
    """
    # Create a new figure with a specified size
    fig = plt.figure(figsize=(16, 6))

    # Create subplots: one for 2D plot and one for 3D plot
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, projection='3d')

    # Plot 2D lines for sudiom and potasiom
    ax1.plot(df['time'], df['sudiom'], label='Sudiom', color='red')
    ax1.plot(df['time'], df['potasiom'], label='Potasiom', color='blue')

    # Set labels and title for the 2D plot
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Values')
    ax1.set_title('Sudiom and Potasiom')

    # Plot 3D points based on the mean values of sudiom and potasiom
    if df['sudiom'].mean() > df['potasiom'].mean():
        ax2.plot(df['time'], df['sudiom'], zs=df['potasiom'], zdir='y', c='r', marker='o', label='Sudiom')
    else:
        ax2.plot(df['time'], df['potasiom'], zs=df['sudiom'], zdir='y', c='b', marker='x', label='Potasiom')

    # Set labels and title for the 3D plot
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Sudiom')
    ax2.set_zlabel('Potasiom')
    ax2.set_title('Brain Data (3D)')

    # Add a legend to the 3D plot
    ax2.legend()

    # Show the plots
    plt.show()

# Generate Brain data and plot it
DF_3D = Brain_Data()
Brain_Plot_3D(DF_3D)








"""
Created on Tue Mar 26 12:04:01 2024

@author: Nasim Heidaran
"""






'''

____________________________Rejection Calculation______________________________

Mathematical modeling of porous (like UF membranes) and dense membranes (like
RO) is a relatively simple procedure. But in the case of NF membrane, modeling 
becomes alittle complicated due to their special structures.
Several mechanisms are suggested to explain the nanofilteration performance of 
these types of membranes.
In this code 3 mechanism are considered to predict rejections of 6 membranes.
the calculated values are compared to experimental ones and the model with the
least average error is introduced as the best model.

 1.Sieving(Ferry-Rankin equation)
 Ferry-Rankin equation is used for calculating membrane rejection.
 The key parameters in this equation are membrane pore radii and solute radii.
 
 2.Hindrance Transport (DSPM-DE equation)
 In this model the role of membrane pore wall is also important in separation.
 The key parameters in this equation are membrane pore radii and solute radii.
 The rejection is calculated by Simplified DSPM-DE equation
 
 3.Donnan Effect (Spiegler_Kedem)
 In this model the key parameter is membrane pore charge density which is
 calculated from zeta potential experimental values.The rejection is 
 calculated by Spiegler_Kedem equation 
 
_______________________________Flux Calculation________________________________

In this code,  Hagen-Poiseuille equation is used for calculating porous 
membrane flux. and Darcy law is used for dense ones.
The key parameters in this equation are membrane pore radii, membrane porosity
thickness of membrane, and dynamic viscosity of solution. In this 
code it is sopposed to calculate volumetric flux for 20 membranes with 
different pore sizes and porosities.  

_____________________Concentration Polarization Evaluation_____________________

Concentration polarization occurs when the concentration of solute increases at
the boundary layer close to the membrane surface. So calculating solute 
concentration on membrane surface is the aim of this code which uses film thory.
'''
#____________________________Rejection Calculation_____________________________

#importing required libraries
    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

#Importing experimental datas
Experimental_Datas=pd.read_excel('C://Users//nasim//Documents//python//Rejection-Flux.xlsx')

#operating condition of filteration
Dox_Size=1.68*10**-9 #m (Solute radius)
Pump_Pressure=4*10**5 #pa (Pressure of pump)
Bulk_Concentration=9.2*10**-5 #mol/m3 (Concentration of solute in polluted water)

#experimental Datas:
Pore_Sizes=np.array(Experimental_Datas['Membrane Pore Size '])  
Experimental_Rejection=np.array(Experimental_Datas['Experimental Rejection'])
Experimental_Flux=np.array(Experimental_Datas['Experimental Flux'])
Experimental_Zeta=np.array(Experimental_Datas['Mean Zeta Potential (mV)'])
Membrane_Code=np.array(Experimental_Datas['Membrane Name'])


def Membrane_Rejection(Ex_R,Ex_J,rs,rp,Zeta,Cm,DelP):
    #constants:
    D=10**-9
    RR=8.314
    T=298
    
    #rejection calculation
    Calculated_Rejection=np.array([(1-2*(1-(rs/rp)**2)+(1-(rs/rp)**4))]) #Ferry-Rankin equation
   
    #error evaluation:
    errorr=abs((Ex_R-Calculated_Rejection)/Ex_R)
    error=np.array(errorr)
    Mean_Error=np.mean(error)
    if Mean_Error<0.15: #if mean error is less than 15% the assumption that the separation performance is due to sieving mechanism, is true
        return Calculated_Rejection
        print('Separation performance is due to sieving mechanism')

    else:
        theta=rs/rp
        mean_theta=np.mean(theta)
        Ka=(1+3.867*theta-1.907*theta**2-0.834*theta**3)/(1+1.867*theta-0.741*theta**2) #convection steric hindreance factor
        
        a=[]
        for i in theta:
            a.append(math.log10(i))
        a1=np.array(a)    
        Kd=((1+(9/8)*a1*theta-1.56*theta+0.53*theta**2+1.95*theta**3-2.82*theta**4+0.27*theta**5+1.1*theta**6-0.44*theta**7)/((1-theta)**2)) ##convection steric hindreance factor
        
        phiSi=(1-(rs/rp))**2
        
        exjka=[]
        JKa=Ex_J*Ka/Kd
        for i in JKa:
            exjka.append(np.exp(-0.11*10**6*i))
      
        #rejection calculation
        Calculated_Rejection2=1-((Ka*phiSi)/(1-(1-Ka*phiSi)*exjka)) #Simplified DSPM-DE Model
        
        #error evaluation:
        errorrr=abs((Ex_R-Calculated_Rejection2)/Ex_R)
        error2=np.array(errorrr)
        Mean_Error2=np.mean(error2)
        if Mean_Error2<0.15:
            print('Separation performance is due to Sieving Mechanism & Hindrance Transport')
            return Calculated_Rejection2
        else:
            Xd=-0.00027+0.001*Zeta
            
            plt.plot(Xd,Zeta,marker='o',ms=10,mec='palevioletred',mfc='palevioletred',color='TEAL',linewidth=5)
            
            ax = plt.gca()
            ax.set_facecolor('ghostwhite')

            font_t={'family':'serif','color':'black','size':15}
            font_x={'family':'serif','color':'black','size':10}
            font_y={'family':'serif','color':'black','size':10}

            plt.title('Zeta Potential vs Xd',fontdict=font_t)
            plt.xlabel('Membrane Pore Density (mol/lit)',fontdict=font_x)
            plt.ylabel('Zeta Potential (mV)',fontdict=font_y)
            
            
            etha=abs(Xd/Cm)
            sigma=1-2/((etha)+(etha**2+4)**0.5)
            omega=D/(RR*T)*(1-sigma)
            Pe=0.14*10**-3*Ex_J/D
           
            FF=[]
            n=(-(1-sigma)*D*Pe/(RR*T*omega))
            for i in n:
                FF.append(np.exp(i))
            F=np.array(FF)
            
            #rejection calculation
            Calculated_Rejection3=sigma*(1-F)/(1-sigma*F) #Spiegler-Kedem equation           
            print('Separation performance is due to Donnan Effect')
            return Calculated_Rejection3

Rejection=Membrane_Rejection(Experimental_Rejection,Experimental_Flux,Dox_Size,Pore_Sizes,Experimental_Zeta,Bulk_Concentration,Pump_Pressure)

print (Rejection)

eror=abs(Experimental_Rejection-Rejection/Experimental_Rejection)
eror=np.array(eror)
Mean_Eror=np.mean(eror)*100

m='average error in this model is '+ str(Mean_Eror) +' percent'
print(m)



#________________________________Flux Calculation______________________________

#Importing experimental datas
Experimental_Datas=pd.read_excel('C://Users//nasim//Documents//python//Flux.xlsx')
Pump_DeltaP=4*10**5
Diffusivity=10**-9
Thickness=0.11*10**-3
Bulk_Concentration=9.2*10**-5
Pore_Size=np.array(Experimental_Datas['Pore Size'])
Porosity=np.array(Experimental_Datas['Porosity'])
Membrane_Code=np.array(Experimental_Datas['Membrane Code']) 
Permeat_Concentration=np.array(Experimental_Datas['Permeate Concentration'])

def Membrane_Flux(r,e,p,l,Cp,Cb):
    mean_pore_size=np.mean(Pore_Size)
    if mean_pore_size<1*10**-9:
        J=(D*(Cb-Cp))/l
        print('Membrane is dense and permeability should be calculated based on Darcy law.')
    else:
        J=(r**2*e*p)/(8*10**-3*l)
        print('Membrane is porous and permeability should be calculated based on Hagen–Poiseuille equation.')
    return J
Flux=Membrane_Flux(Pore_Size,Porosity,Pump_DeltaP,Thickness,Permeat_Concentration,Bulk_Concentration)

plt.bar(Membrane_Code,Flux,color='mediumvioletred', width=0.3)


font_ttt={'family':'serif','color':'black','size':15}
font_xxx={'family':'serif','color':'black','size':10}
font_yyy={'family':'serif','color':'black','size':10}

plt.title('Mmebrane Flux',fontdict=font_ttt,pad=25)
plt.xlabel('Membrane code',fontdict=font_xxx)
plt.ylabel('Mmebrane Pure Water Flux (m/s) ',fontdict=font_yyy)
plt.show()
#________________________________Film Treory___________________________________

#Importing experimental datas
Experimental_Datas=pd.read_excel('C://Users//nasim//Documents//python//Film Theory.xlsx')

Membrane_length=0.09
Mmembrane_width=0.04 
Mmebrane_thickness=0.11 *10**-3
Velocity=0.5
Bulk_Concentration=9.2*10**-5
Permeate_Concentration=np.array(Experimental_Datas['Permeate Concentration']) 
Flux=np.array(Experimental_Datas['Flux'])
Membrane_Code=np.array(Experimental_Datas['Membrane Code']) 

def Membrane_Concentration(a,b,u,J,Cb,Cp,l):
    dh=4*a*b/(a+b)
    Re=1000*u*dh/(10**-3)
    Sc=10**-3/(1000*10**-9)
    k=((0.664*Re**0.5*Sc**0.33*(dh/l)**0.5)*10**-9)/dh
   
    v=[]
    Jk=J/k
    for i in Jk:
        v.append(np.exp(i))
    
    Cm=(Cb-Cp)*v+Cp
    Polarization=(Cm-Cb)/Cb
    mean_Polarization=np.mean(Polarization)
    
    if mean_Polarization<0.05:
        print('Concentration polarization is neglible. So, membrane has a low probability for fouling.')
    else:
        print('Concentration polarization is high. So, fouling may occures.')
    

    return(Cm)
    
Membrane_Surface_Concentration=Membrane_Concentration(Membrane_length,Mmembrane_width,Velocity,Flux,Bulk_Concentration,Permeate_Concentration,Mmebrane_thickness)

plt.bar(Membrane_Code, Membrane_Surface_Concentration,color='mediumvioletred', width=0.3)


font_tt={'family':'serif','color':'black','size':15}
font_xx={'family':'serif','color':'black','size':10}
font_yy={'family':'serif','color':'black','size':10}

plt.title('Mmebrane Surface Concentration',fontdict=font_tt,pad=20)
plt.xlabel('Membrane code',fontdict=font_xx)
plt.ylabel('Cm',fontdict=font_yy)
plt.show()




import numpy as np

import matplotlib.pyplot as plt

import pandas as pd

import math

from matplotlib.patches import Ellipse, Polygon



f_loc='C://Users//Haj Abedi//Dropbox//Nasim//Python//Project//Data.xlsx'

def SI_Calculation(P,PC,Density=1):

    

    '''

    This function is used for Separation Index Calculation

    P : Pressure (bar)

    Density : Feed Density(g/cm3)

    PC :  Pollutant concentration in Feed (g/L)

    Returns Separation Index and Flux & Rejection & Rejection Charts

    '''

    

    Data=pd.read_excel(f_loc)

    Data.columns

    J=Data['Flux']

    R=Data['Rejection']



    SI=(Density-(1-R)*PC/1000)*(J/(P*((1-R)**0.41)))

    Mem_Code=np.array(Data['Mem Code'])

    Flux=np.array(Data['Flux'])

    Rejection=np.array(Data['Rejection'])



    font={'family':'serif','color':'k','size':'20'}

    

    c=np.array([])

    for i in range (0,len(Flux)):

        if Flux[i]<100:

            a=np.array(['.'])

            c=np.concatenate((c,a))

        elif Flux[i]<200:

            a=np.array(['o'])

            c=np.concatenate((c,a))

        else:

            a=np.array(['O'])

            c=np.concatenate((c,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Flux','o:Medium Flux','O:High Flux']

    for i in range(0,len(Flux)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(Mem_Code,Flux,color='w',edgecolor='c',hatch=c,linewidth=1,yerr=10,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Flux Chart',fontdict=font)

    plt.xlabel('Membrane code',fontdict=font)

    plt.ylabel('Flux',fontdict=font)

    ax.legend(title='Flux Range')

    plt.show()

    

    d=np.array([])

    for i in range (0,len(Rejection)):

        if Rejection[i]<0.6:

            a=np.array(['.'])

            d=np.concatenate((d,a))

        elif Rejection[i]<0.75:

            a=np.array(['o'])

            d=np.concatenate((d,a))

        else:

            a=np.array(['O'])

            d=np.concatenate((d,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Rejection','o:Medium Rejection','O:High Rejection']

    for i in range(0,len(Rejection)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(Mem_Code,Rejection,color='w',edgecolor='c',hatch=d,linewidth=1,yerr=0.01,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Rejection Chart',fontdict=font)

    plt.xlabel('Membrane code',fontdict=font)

    plt.ylabel('Rejection',fontdict=font)

    ax.legend(title='Rejection Range')

    plt.show()

    

    f=np.array([])

    for i in range (0,len(SI)):

        if SI[i]<250:

            a=np.array(['.'])

            f=np.concatenate((f,a))

        elif SI[i]<500:

            a=np.array(['o'])

            f=np.concatenate((f,a))

        else:

            a=np.array(['O'])

            f=np.concatenate((f,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low SI','o:Medium SI','O:High SI']

    for i in range(0,len(SI)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(Mem_Code,SI,color='w',edgecolor='c',hatch=f,linewidth=1,yerr=10,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('SI Chart',fontdict=font)

    plt.xlabel('Membrane code',fontdict=font)

    plt.ylabel('SI',fontdict=font)

    ax.legend(title='SI Range')

    plt.show()    

    

    return SI

SI=SI_Calculation(0.6,50/1000)





e_loc='C://Users//Haj Abedi//Dropbox//Nasim//Python//Project//Porosity Data.xlsx'

def Porosity(Density=1):

    

    '''

    Ww : weight of wet samples (g)

    Wd : weight of dry samples (g)

    V : denotes the sample volume (cm3)

    Density : is the water density (g/cm3).The default is 1.

    Returns the porosity of membranes

    '''

    

    Porosity_Data=pd.read_excel(e_loc)

    Porosity_Data.columns

    Ww=Porosity_Data['Ww']

    Wd=Porosity_Data['Wd']

    V=Porosity_Data['V']

    

    Porosity=(Ww-Wd)/(Density*V)

    membrane=np.array(Porosity_Data['membrane'])



    font={'family':'serif','color':'k','size':'20'}

    

    c=np.array([])

    for i in range (0,len(Porosity)):

        if Porosity[i]<0.9:

            a=np.array(['.'])

            c=np.concatenate((c,a))

        elif Porosity[i]<1:

            a=np.array(['o'])

            c=np.concatenate((c,a))

        else:

            a=np.array(['O'])

            c=np.concatenate((c,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Porosity','o:Medium Porosity','O:High Porosity']

    for i in range(0,len(Porosity)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(membrane,Porosity,color='w',edgecolor='c',hatch=c,linewidth=1,yerr=0.05,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Porosity Chart',fontdict=font)

    plt.xlabel('membrane',fontdict=font)

    plt.ylabel('Porosity',fontdict=font)

    ax.legend(title='Porosity Range')

    plt.show()

    

    return Porosity

Porosity=Porosity()



def Tortuosity():

    

    '''

    Returns the Pore Tortuosity of membranes

    '''

    

    Porosity_Data=pd.read_excel(e_loc)

    Tortuosity=((2-Porosity)**2)/Porosity

    membrane=np.array(Porosity_Data['membrane'])



    font={'family':'serif','color':'k','size':'20'}

    

    c=np.array([])

    for i in range (0,len(Tortuosity)):

        if Tortuosity[i]<0.75:

            a=np.array(['.'])

            c=np.concatenate((c,a))

        elif Tortuosity[i]<1.25:

            a=np.array(['o'])

            c=np.concatenate((c,a))

        else:

            a=np.array(['O'])

            c=np.concatenate((c,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Tortuosity','o:Medium Tortuosity','O:High Tortuosity']

    for i in range(0,len(Tortuosity)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(membrane,Tortuosity,color='w',edgecolor='c',hatch=c,linewidth=1,yerr=0.05,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Tortuosity Chart',fontdict=font)

    plt.xlabel('membrane',fontdict=font)

    plt.ylabel('Tortuosity',fontdict=font)

    ax.legend(title='Tortuosity Range')

    plt.show()

    

    return Tortuosity

Tortuosity=Tortuosity()





g_loc='C://Users//Haj Abedi//Dropbox//Nasim//Python//Project//Pore Size Data.xlsx'

def Pore_Size(A,P,Vis=8.9*1e-4):

    

    '''

    A=shows the membrane effective surface area (m2)

    P : indicates the utilized operational pressure (Pa)

    Vis : represents the water viscosity (Pa⋅s). The default is 8.9*1e-4.

    Returns the Pore Size of membranes in nm

    '''

    

    Pore_Size_Data=pd.read_excel(g_loc)

    Pore_Size_Data.columns

    q=Pore_Size_Data['q']

    l=Pore_Size_Data['l']



    Pore_Size=((2.9-1.75*Porosity)*(8*Vis*q*l/1000)/(Porosity*A*P))**(0.5)*1E9

    membrane=np.array(Pore_Size_Data['membrane'])



    font={'family':'serif','color':'k','size':'20'}

    

    c=np.array([])

    for i in range (0,len(Pore_Size)):

        if Pore_Size[i]<4.5:

            a=np.array(['.'])

            c=np.concatenate((c,a))

        elif Pore_Size[i]<5.5:

            a=np.array(['o'])

            c=np.concatenate((c,a))

        else:

            a=np.array(['O'])

            c=np.concatenate((c,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Pore_Size','o:Medium Pore_Size','O:High Pore_Size']

    for i in range(0,len(Pore_Size)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    ax.bar(membrane,Pore_Size,color='w',edgecolor='c',hatch=c,linewidth=1,yerr=0.05,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Pore_Size Chart',fontdict=font)

    plt.xlabel('membrane',fontdict=font)

    plt.ylabel('Pore Size',fontdict=font)

    ax.legend(title='Pore_Size Range')

    plt.show()

    

    return Pore_Size

Pore_Size=Pore_Size(0.0024,400000)





def Pore_Number(A):

    

    '''

    A=shows the membrane effective surface area (m2)

    Returns the Pore Number of membranes

    '''

    

    Pore_Number=(Porosity*A)/((math.pi)*(((Pore_Size*1E-9)**2)/4))



    font={'family':'serif','color':'k','size':'20'}

    

    c=np.array([])

    for i in range (0,len(Pore_Number)):

        if Pore_Number[i]<1*1E14:

            a=np.array(['.'])

            c=np.concatenate((c,a))

        elif Pore_Number[i]<2*1E14:

            a=np.array(['o'])

            c=np.concatenate((c,a))

        else:

            a=np.array(['O'])

            c=np.concatenate((c,a))

    fig, ax = plt.subplots()

    

    bar_labels=['.:Low Pore_Number','o:Medium Pore_Number','O:High Pore_Number']

    for i in range(0,len(Pore_Number)-3):

        m=['_.']

        bar_labels=bar_labels+m

        

    Pore_Size_Data=pd.read_excel(g_loc)   

    membrane=np.array(Pore_Size_Data['membrane'])

    ax.bar(membrane,Pore_Number,color='w',edgecolor='c',hatch=c,linewidth=1,yerr=0.05,ecolor='c',width=0.85,label=bar_labels)   

    plt.title('Pore_Number Chart',fontdict=font)

    plt.xlabel('membrane',fontdict=font)

    plt.ylabel('Pore Number',fontdict=font)

    ax.legend(title='Pore_Number Range')

    plt.show()

    

    return Pore_Number

Pore_Number=Pore_Number(0.0024)

gfhfg=1



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''
This function is used for calculating rejection and flux of nanofilteration 
membrane for separation of anti-tumor drug Doxorubicin.z

_____________Rejection Calculation_____________

In this code, Ferry-Rankin equation is used for calculating membrane rejection.
The key parameters in this equation are membrane pore radii and solute pore   
radii. In this code it is sopposed to calculate rejection for 22 membranes with 
different pore sizes.

_____________Flux Calculation_____________

In this code,  Hagen-Poiseuille equation is used for calculating membrane flux.
The key parameters in this equation are membrane pore radii, membrane pore 
numbers, thickness of membrane, and dynamic viscosity of solution. In this 
code it is sopposed to calculate volumetric flux for 22 membranes with 
different pore sizes and pore numbers.  
'''
#_____________________________________________________________________________
#importing required libraries:
    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import solve_bvp
import scipy.optimize
import math


Experimental_Datas=pd.read_excel('C://Users//nasim//Documents//python//Rejection-Flux.xlsx')
Dox_Size=1.68*10**-9 #m

def Membrane_Rejection(Experimental_Datas):
    Pore_Sizes=np.array(Experimental_Datas['Membrane Pore Size ']) 
    Calculated_Rejection=(1-2*(1-(Dox_Size/Pore_Sizes)**2)+(1-(Dox_Size/Pore_Sizes)**4))
    Experimental_Rejection=np.array(Experimental_Datas['Experimental Rejection'])
    error=abs((Experimental_Rejection-Calculated_Rejection)/Experimental_Rejection)
    Mean_Error=np.mean(error)
    if Mean_Error<0.05:
        Rejection=Calculated_Rejection
        return(Rejection)
        print('Separation performance is due to sieving mechanism')
    else:
        Zeta_Potential=np.array(Experimental_Datas['Mean Zeta Potential (mV)'])
        Xd=-0.27+2.8*Zeta_Potential #a corralation between zeta potential and membrane charge density. it is worth noty that zeta potential is a measurable parameter while,membrane charge density is not.
        theta = (1-Dox_Size/Pore_Sizes)**2
        Kc =(2-theta)*(1+0.054*Dox_Size/Pore_Sizes-0.988*(Dox_Size/Pore_Sizes)**2+0.441*(Dox_Size/Pore_Sizes)**3)
        Kd=1-2.3*Dox_Size/Pore_Sizes+1.15*(Dox_Size/Pore_Sizes)**2+0.224*(Dox_Size/Pore_Sizes)**3
        C0=9.2*10**-5
        Cp=C0
       
        def ODE_Function(x,c,Cp): #Define simultaneous ODE equations
            Experimental_J=np.array(Experimental_Datas['Experimental Flux']) 
            z = 1
            D=10**-9
            Dp=Kd*D
            dcdx = Experimental_J/Dp*(Kc*c-Cp)-z*c*(Experimental_J*z/Dp*(Kc*c-Cp))/((z**2*c))
            return dcdx
            
        def Boundary_Conditions(ca,cb,Cp):
            T=298
            F=96487
            R=8.314
            z=1
        
            a1=(abs(Xd)/(z*C0*theta))
            phiDa=[]
            for i in a1:
                phiDa.append(math.log10(i)/(F*z/R/T))
            phiDa_array=np.array(phiDa)
            c0_m=C0*theta*np.exp(-F*z*phiDa_array/R/T)
            
            a2=(abs(Xd)/(z*Cp*theta))
            phiDb=[]
            for i in a2:
                phiDb.append(math.log10(i)/(F*z/R/T))
            phiDb_array=np.array(phiDb)
            cb_m=Cp*theta*np.exp(-F*z*phiDb_array/R/T)
                        
            return np.array([c0_m - ca[0], cb_m - cb[-1]])
        x = np.linspace(0, 0.14,6)
        y= np.zeros((1, x.size))
        p = np.array([Cp])
        res= solve_bvp(ODE_Function,Boundary_Conditions,x,y,p=p)
        x_plot = np.linspace(0, 0.14, 100)
        res_plot=res.sol(x_plot)[1]
        cp= res_plot[-1]
        Rejection=1-cp/C0
        return(Rejection)
Rejection=Membrane_Rejection(Experimental_Datas)

'''
yeki az mabahese mortabeT ba Industrial Engineering mabhase Amare tosifi mibashad
amare tosifi shamele:
    1. amare haye tosif konande marboot be goroohi az dade ha shamele:
        1.1. amare haye tamrkozi mandande mean, median, mod
        1.2. amare haye parakandegi hamchon standard deviation, variance
        1.3. amare haye shekle tozi hamchon skewness
    2. nemoodar haye tosifi hamchon:
        2.1. Scatterplot
        2.2. Histogram
        2.3. Boxplot
        and so on ...
        
dar in barname 2 tabe tarif mikonim:
    
    1.Des_Stat: in tabe ettela,ate marboot be amare tosifi ra baraye yek series moshakhas erae midahad

    2.Des_Plot: in tabe nemoodare morede nazare user marboot be 2 series ra rasm mikondad
'''


Finance_Data=pd.read_excel('C://Users//user//Desktop//Financial.xlsx')
#calling data sample via Pandas library

Random_Data=pd.DataFrame(np.array([[1,2,3],[4,5,6],[7,8,9]]),columns=['A','B','C'])
#calling data in a reqular way!

MSP_Class=pd.DataFrame(
    {
     'Name':['Mojtaba Bakhshi','Amir Sabri','Mohammad Bakhtiyari'],
     'Age':[25,24,23],
     'Sex':['Male','Male','Male'],
     'Education':['Bs','Ms','Bs']
     }
    )
#calling data using dictionaries. When using a Python dictionary of lists, the dictionary keys will be used as column headers and the values in each list as columns of the DataFrame.

def Des_Stat ():
    
    Column=Finance_Data[input('please tell me which column you interested it: ')]
    Descriptive_Statistical=Column.describe()
    
    return Descriptive_Statistical

#for example

Des_Stat()

'''
Out[10]: 
count    700.000000
mean     118.428571
std      136.775515
min        7.000000
25%       12.000000
50%       20.000000
75%      300.000000
max      350.000000
Name: Sale Price, dtype: float64
'''
 # in tabe name sotooni ra ke qasde daryafte ette'lat ra darim ra migirad va az dakhele dataframe ettela,ate marboot be an sotoon ra dakhele variable Column zakhire mikonad
 # function 'Describe()' yeki az tabe haye ketabkhane pandas mibashad ke kollie ette'laate marboot be amare tosifi ra ejra mikonad
 
 
def Des_Plot () :
    
    Xaxis=Finance_Data[input('Enter the name of x_axis column: ')]
    Yaxis=Finance_Data[input('Enter the name of y_axis column: ')]
    plot=input('Which plot would you prefer to show: ')
    
    if plot=='Scatter':
        plt.scatter(Xaxis, Yaxis)
        plt.title(plot, c='purple', fontsize=24, fontname='Times New Roman')
        plt.xlabel(Xaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.ylabel(Yaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.show()
    elif plot=='Histogram':
        plt.hist2d(Xaxis, Yaxis)
        plt.title(plot, c='purple', fontsize=24, fontname='Times New Roman')
        plt.xlabel(Xaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.ylabel(Yaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.show()
    elif plot=='Plot':
        plt.plot(Xaxis, Yaxis)
        plt.title(plot, c='purple', fontsize=24, fontname='Times New Roman')
        plt.xlabel(Xaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.ylabel(Yaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.show()
    elif plot=='Stem':
        plt.stem(Xaxis, Yaxis)
        plt.title(plot, c='purple', fontsize=24, fontname='Times New Roman')
        plt.xlabel(Xaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.ylabel(Yaxis.name, c='crimson', fontsize=14, fontname='Times New Roman')
        plt.show()
    else:
        print('Maybe you mean a plot that I cant show you!')

Des_Plot()










import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def XRD_Analysis(file,which,peak=0):
    '''
    

    Parameters
    ----------
    file : str
        the variable in which you saved the .cvs file path         
    which : str
        which operation you want to perform on the file      
    peak : float, optional
        2θ for the peak you want to analyse. The default is 0.     

    Returns
    -------
    fwhm : float
        value of FWHM for the peak you specified.

    '''
    
    df=pd.read_csv(file)
    npar=pd.DataFrame.to_numpy(df)

    if which=='plot':
        angle=df['angle']
        intensity=df['intensity']
        plt.plot(angle,intensity,color='k')
        font_title={'family':'serif','color':'blue','size':20}
        plt.title('XRD pattern',fontdict=font_title)
        font_label={'family':'times new roman','color':'black','size':15}
        plt.xlabel('angle (2θ)',fontdict=font_label)
        plt.ylabel('intensity (a.u.)',fontdict=font_label)
        plt.grid(axis='x',which='both')
        plt.xticks(np.arange(0,max(angle),5))
        plt.xlim([np.min(npar,axis=0)[0], np.max(npar,axis=0)[0]])
        plt.yticks([])
        plt.ylim([0, 1.1*np.max(npar,axis=0)[1]])
        plt.tick_params(axis='x',direction='in')
        plt.show()
        return None
    elif which=='fwhm':
        diff=int((npar[1,0]-npar[0,0])*1000)/2000
        for i in range(int(len(npar)/2)+1):
            if -diff<npar[i,0]-peak<diff:
                pl=i
                ph=i
                p=i
                break
        while pl>0:
            if ((npar[pl,1]-npar[pl-1,1])/(npar[pl-1,1]-npar[pl-2,1]))>1.04 and (npar[pl-1,1]-np.min(npar,axis=0)[1])/(np.max(npar,axis=0)[1]-np.min(npar,axis=0)[1])<0.4:
                in_low_1=npar[pl-1,1]
                break
            pl=pl-1
        while ph>0:
            if ((npar[ph+2,1]-npar[ph+1,1])/(npar[ph+1,1]-npar[ph,1]))<0.96 and (npar[ph+1,1]-np.min(npar,axis=0)[1])/(np.max(npar,axis=0)[1]-np.min(npar,axis=0)[1])<0.4:
                in_low_2=npar[ph+1,1]
                break
            ph=ph+1
        in_low=(in_low_1+in_low_2)/2
        h=npar[p,1]-in_low
        hm=in_low+h/2
        diff_in=[]
        hm_i=[]
        for l in range(len(npar)-1):
            diff_in.append((npar[l+1,1]-npar[l,1])/2)
        for j in range(2):
            for k in range(int(len(npar)/2)+1):
                c=((-1)**j)*k
                if abs(npar[p+c,1]-hm)<abs(max(diff_in)):
                    hm_i.append(p+c)
                    break
        fwhm=npar[hm_i[0],0]-npar[hm_i[1],0]
        return fwhm
    else:
        print('The which argument not valid')
        return None



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Loading data
mydf=pd.read_csv('nanofiber-stress-strain.csv')
mydf1=mydf[39:175]
mydf2=mydf[203:414]
mydf3=mydf[442:]

# Defining a function to change column names and clean data
def preprocess(df):
    # Rename columns
    df.columns = df.iloc[0]
    df = df[1:]
    df.rename(columns={'Tensile stress MPa': 'Stress (MPa)', 'Tensile strain %': 'Strain (%)'}, inplace=True)
    
    # Reseting index
    df.reset_index(inplace=True, drop=True)
   
    #changing the type of the columns
    df = df.astype(float)
    return df

# applying preprocess function on dataframes
df_new1 = preprocess(mydf1)
df_new2 = preprocess(mydf2)
df_new3 = preprocess(mydf3)

def Stress_Strain_Curve(input_data, action):
    stress = input_data['Stress (MPa)']
    strain = input_data['Strain (%)']
    
    if action == 'plot':
        # Plotting data
        fig, ax = plt.subplots(figsize=(8, 6))
        plt.plot(strain, stress, linewidth=2, color='royalblue', marker='o', markersize=5, label='Stress-Strain Curve')
        plt.title('Stress-Strain Curve', fontsize=16)
        plt.xlabel('Strain (%)', fontsize=14)
        plt.ylabel('Stress (MPa)', fontsize=14)
        plt.xlim([0, strain.max()])
        plt.ylim([0, stress.max()])
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()
   
    elif action == 'max stress':
        # Calculation of the maximum stress
        stress_max = stress.max()
        return stress_max
    
    elif action == 'young modulus':
        # Calculation of Young's Modulus
        slope_intercept = np.polyfit(strain, stress, 1)
        return slope_intercept[0]
    
# Example:
stress_max = Stress_Strain_Curve(df_new1, 'max stress')
young_modulus = Stress_Strain_Curve(df_new1, 'young modulus')
Stress_Strain_Curve(df_new1, 'plot')
    
print('Maximum Stress (MPa):', stress_max)
print("Young's Modulus (MPa):", young_modulus)



import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

def aerospace (CSV,which):
    '''
    this function has the ability to convert your datas into 
    answers that you need 
    your datas should be in Newton and M**2 format 
    in this program we will be using presure as the output data 
    if you want to make a sketch Use === Plot
    if you want to check the max_presure use === MaxPer
    '''
    mydf = pd.read_csv(CSV)
    mydff = np.array(mydf)
    mydf1 = pd.DataFrame(mydff,columns=['Newton','Area'])
    mydf2 = mydf1['Newton']/mydf1['Area']
    mydf3 = pd.concat(mydf1,mydf2)
    if which == 'Plot':
        plt.plot(mydf1['Newton'],mydf1['Area'])
        plt.xlabel('Area')
        plt.ylabel('Newton')
        plt.show()
    if which == 'MaxPer':
        max_p = mydf3.max()
        return max_p
        

"""
Created on Sun Sep  2 12:11:53 2023

@author: ALirezaPeymani


Remove spikes from Raman spectra of the Polybutadiene by the Z-scores  algorithm 

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# load the data as data frame
df = pd.read_csv("C://Users//ALPHA65//Quiz//Ram_PBDEN.csv")
# Transform the data to a numpy array
wavelength = df['Wavelength']
intensity = df['Intensity']
# Plot the spectrum:
plt.plot(wavelength, intensity)
plt.title('Spectrum', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('Wavelength (cm-1)', fontsize = 20)
plt.ylabel('Intensity (a.u.)', fontsize = 20)
plt.show()


#calculate the z-scores for the points in the Raman spectrum
#z-scores = (x(i)-μ) / σ
#where μ is the mean and σ is the standard deviation of the population x 
#(x(i) represent the values of a single Raman spectrum)
#The z-scores tell how far a value is from the average in units of standard deviation.

def z_score(intensity):
 mean_int = np.mean(intensity)
 std_int = np.std(intensity)
 z_scores = (intensity - mean_int) / std_int
 return z_scores 

intensity_z_score = np.array(z_score(intensity))
plt.plot(wavelength, intensity_z_score)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel( 'Wavelength (cm-1)' ,fontsize = 20)
plt.ylabel( 'Z-Score' ,fontsize = 20)
plt.show()   


def z_score(intensity):
 mean_int = np.mean(intensity)
 std_int = np.std(intensity)
 z_scores = (intensity - mean_int) / std_int
 return z_scores
threshold = 3.5
intensity_z_score = np.array(abs(z_score(intensity)))
plt.plot(wavelength, intensity_z_score)
plt.plot(wavelength, threshold*np.ones(len(wavelength)), label = 'threshold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel( 'Wavelength' ,fontsize = 20)
plt.ylabel( '|Z-Score|' ,fontsize = 20)
plt.show()


threshold = 3.5
# 1 is assigned to spikes, 0 to non-spikes:
spikes = abs(np.array(z_score(intensity))) > threshold
plt.plot(wavelength, spikes, color = 'red')
plt.title('Spikes:'  + str(np.sum(spikes)), fontsize = 20)
plt.grid()
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel( 'Wavelength' ,fontsize = 20)
plt.ylabel( 'Z-scores > ' + str(threshold) ,fontsize = 20)
plt.show()


#The second approach  modified Z-score method uses the median (M) and median absolute deviation (MAD) rather than the mean and standard deviation:

#z(i) = 0.6745 (x(i)-M) / MAD
#The multiplier 0.6745 is the 0.75th quartile of the standard normal distribution.


def modified_z_score(intensity):
 median_int = np.median(intensity)
 mad_int = np.median([np.abs(intensity - median_int)])
 modified_z_scores = 0.6745 * (intensity - median_int) / mad_int
 return modified_z_scores
threshold = 3.5
intensity_modified_z_score = np.array(abs(modified_z_score(intensity)))
plt.plot(wavelength, intensity_modified_z_score)
plt.plot(wavelength, threshold*np.ones(len(wavelength)), label = 'threshold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('Wavelength (cm-1)' ,fontsize = 20)
plt.ylabel('|Modified Z-Score|' ,fontsize = 20)
plt.show()



# 1 is assigned to spikes, 0 to non-spikes:
spikes = abs(np.array(modified_z_score(intensity))) > threshold
plt.plot(wavelength, spikes, color = 'red')
plt.title('Spikes: ' + str(np.sum(spikes)), fontsize = 20)
plt.grid()
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel( 'Wavelength' ,fontsize = 20)
plt.ylabel( 'Modified Z-scores > ' + str(threshold) ,fontsize = 20)
plt.show()


#The third approach
#The third approach propose to make advantage of the high intensity and small width of spikes.  
#،herefore use the difference between consecutive spectrum points Dx(i) = x(i)-x(i-1) to calculate the z-scores.

# First we calculated Dx(i):
dist = 0
delta_intensity = [] 
for i in np.arange(len(intensity)-1):
 dist = intensity[i+1] - intensity[i]
 delta_intensity.append(dist)
delta_int = np.array(delta_intensity)
# Alternatively to the for loop one can use: 
# delta_int = np.diff(intensity)

intensity_modified_z_score = np.array(modified_z_score(delta_int))
plt.plot(wavelength[1:], intensity_modified_z_score)
plt.title('Modified Z-Score using Dx(i)')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('Wavelength (cm-1)', fontsize = 20)
plt.ylabel('Score', fontsize = 20)
plt.show()

#In order to apply a threshold to exclude spikes, the absolute value of the modified Z-score must be taken:
#|z(i)| =|0.6745 (Dx(i)-M) / MAD|    

threshold = 3.5
intensity_modified_z_score = np.array(np.abs(modified_z_score(delta_int)))
plt.plot(wavelength[1:], intensity_modified_z_score)
plt.plot(wavelength[1:], threshold*np.ones(len(wavelength[1:])), label = 'threshold')
plt.title('Modified Z-Score of Dx(i)', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('Wavelength (cm-1)', fontsize = 20)
plt.ylabel('Score', fontsize = 20)
plt.show()


# 1 is assigned to spikes, 0 to non-spikes:
spikes = abs(np.array(modified_z_score(intensity))) > threshold
plt.plot(wavelength, spikes, color = 'red')
plt.title('Spikes: ' + str(np.sum(spikes)), fontsize = 20)
plt.grid()
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel( 'Wavelength' ,fontsize = 20)
plt.show()


plt.plot(wavelength[1:],
np.array(abs(modified_z_score(delta_int))), color='black', label = '|Modified Z-Score using Dx(i)|')
plt.plot(wavelength, np.array(abs(modified_z_score(intensity))), label = '|Modified Z-Score|', color = 'red')
plt.plot(wavelength, np.array(abs(z_score(intensity))), label = '|Z-Score|', color = 'blue')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel( 'Wavelength(cm-1)' ,fontsize = 20)
plt.ylabel( 'Score' ,fontsize = 20)
plt.legend()
plt.show()




def my_Stress_Strain(input_file,which,count):
    '''
    This function claculates the stress and strain
    Parameters from load and elongation data
    ----------
    input_file : .csv format
        the file must be inserted in csv.
    whcih : str
        please say which work we do ( plot or calculate?).
    count: int
        please enter the yarn count in Tex
    remember: gauge length has been set in 250 mm
    '''

    #convert the file
    mydf=pd.read_csv(input_file)

    if which=='plot':
       
        stress=mydf['Load']/count
        strain=mydf['Extension']/250
        plt.plot(stress,strain)
        plt.title('stress-strain curve')
        plt.xlabel('stress')
        plt.ylabel('strain')
        plt.show()
    
    
    if which=='max stress':
        stress_max=mydf['stress'].max()
        return stress_max

    if which=='max strain':
        strain_max=mydf['strain'].max()
        return strain_max




import pandas as pd
import statistics as st
import numpy as np
import matplotlib.pyplot as plt




def pressure_volume_ideal_gases(file,which):
    '''
    By using this function, the relationship between pressure and volume 
    of ideal gases in thermodynamic will be shown.
    
    input_file : .csv format
        *the file must be inserted in csv.
        
    whhch: str
        what do you want this function to do?
        
    '''
    
    mydf=pd.read_csv(file)
    
    if which=='plot':
        pressure=mydf['pressure']
        volume=mydf['volume']
        
        plt.plot(volume,pressure)
        plt.title('volume_pressure_chart')
        plt.xlabel('pressure')
        plt.ylabel('volume')
        font_title={'family':'serif','color':'black','size':18}
        font_label={'family':'serif','color':'red','size':12}
        plt.show()
        
        
        
    elif which=='min pressure':
        min_pressure=mydf['pressure'].min()
        return min_pressure
    
    elif which=='max pressure':
        max_pressure=mydf['pressure'].max()
        return max_pressure
    
    elif which=='min volume':
        min_volume=mydf['volume'].min()
        return min_volume
    
    elif which=='max volume':
         max_volume=mydf['volume'].max()
         return max_volume
    
    elif which=='average pressure':
         average_pressure=st.mean(pressure)
         return average_pressure
        
    elif which=='average volume':
         average_volume=st.mean(volume)
         return average_volume
            
    elif which=='temperature':
        n=1
        R=0.821
        temperature=(pressure*volume)/(n*R)
        '''
        This formula is from 'Introduction To The thermodinamics Of Materials
        David R. Gaskell'
        
        '''
        return temperature
    
    else:
        print('No answer found')
    
        
    

def Gradient_descent(dataset , which, alpha=0.001):
    '''
    

    Parameters
    ----------
    dataset : str
        path of your data set .
    which : str
        plot or cost(calculuates how much error is there in the prediction).
    alpha : float, optional
        the learning rate(the size of the baby step that we will take for each of the variables). The default is 0.001.


    '''
    data=pd.read_csv(dataset)
    data=data.ffill()
   
    
    if which=='plot':
        plt.xlabel(data.columns[0])
        plt.ylabel(data.columns[1])
        rename_dict={data.columns[0]:'x',data.columns[1]: 'y'}
        data=data.rename(columns=rename_dict)
        x=list(data['x'])
        y=list(data['y'])
        plt.plot(x, y,'o')
        plt.show()
   
    
   
    if which=='cost':
        def cost_function(l,c,p,x,y):
            
            '''
          The linear equation for a line:  y=p*x+c

            Parameters
            ----------
            l : int
                is the number of records.
         x,y,p,c: int
         
           Returns
           -------
         finds the square of the difference..
            '''
            return 1/2/l * sum([(c + p* np.asarray([x[i]]) - y[i])**2 for i in range(l)])
      
        
      
        
      
        c=0
        l = len(data)
        p=1
        oldcost=0
#oldcost is a variable that we will use to save the cost (error) in the last iteration so that we can compare with the current iteration. 
        newcost=0
        gr0=0
#the amount by which we will change the 'c' variable        
        gr1=0
#the amount by which we will change the 'x' variable      
        rename_dict={data.columns[0]:'x',data.columns[1]: 'y'}
        data=data.rename(columns=rename_dict)
        x=data['x']
        y=data['y']
        for i in range(100000000000):
            for j in range(l):
# this 'for loop' essentially finds the average value by which we will change the gradients.                    
                gr0=gr0+((c + p * x[i]) - y[i])
                gr1=gr1+(((c+p* x[i]) - y[i]) * x[i])
                gr0 = gr0/l
                gr1 = gr1/l
# changing the value of 'p' and 'c' to see how it impacts to our cost(error)
                c=c-alpha*gr0
                p=p-alpha*gr1 
              
              
                newcost=cost_function(len(data), c, p, x, y) 
#If the change in cost (error), is less than particular number(here we use 0.001).it has probably reduce the cost(error) as much as it can.
                if abs(oldcost-newcost)<0.001:
                   
                   return newcost
                   break
                else:
                   oldcost=newcost
                   
                   
                   
#ddd=Gradient_descent("C:\\Users\\Yasamin\\Desktop\\pga.csv",'cost')
#df=Gradient_descent("C:\\Users\\Yasamin\\Desktop\\Data1.csv", 'plot')
                   
                   
                    
    
    
    
    

import pandas as pd  


import matplotlib.pyplot as plt   
 
   
def Energie(input_file,which):
    '''
    This is a function to drawing a plot or to calculating 
    the amount of Energie of a Motor to (open/close) a Valve in a cycle, which takes 2.7 secound to open and to close.
    
    ----------
    input_file : .xlsx format
        the file must be inserted in xlsx.
    which : int
        1 : Drawing a Plot
        2 : Calculate the Consupmtion Energie in [mWs]
        please say which work we do ( 1 or 2).

    '''
    
    input_file = 'C:\\Users\\Fust\\Desktop\\Book1.xlsx' 
    
    mydf=pd.read_excel(input_file)
    
    
    if which==1:
        #get the data on each columns
        
        A=mydf['Angle[°]']
        
        Energie =mydf['Energie']
       
        #plotting data
        plt.plot(A, Energie,color = 'green')
        plt.title('Energie of OTT Motor 185000 Cycle')
        plt.xlabel('Angle[°]')
        plt.ylabel('Consupmtion Energie')
        plt.show()
    
    
    if which==2:
        mydf=pd.DataFrame(mydf,columns=['Angle[°]','Power[mW]','Time for a Cycle','Energie'])
        
        summ = mydf['Energie'].sum()                          # The amount of Energie for a half Cycle of Duty in mWs
       
        summ =( summ * 2)/1000                                # The amount of Consumption Energie for a Dutycycle in Ws
        
        return summ
    
Energie( 'C:\\Users\\Fust\\Desktop\\Book1.xlsx', 1)
Energie( 'C:\\Users\\Fust\\Desktop\\Book1.xlsx', 2)




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

adrs='C:\\Users\\Brooz\\Desktop\\datafatigue.csv'
data_fatigue=pd.read_csv('C:\\Users\\Brooz\\Desktop\\datafatigue.csv')
data_fatigue.columns

LoadRange=[]
LoadCycles=[]

for i in range (0,len(data_fatigue)):
    a=data_fatigue['LoadRange LoadCycles'][i].split()
    LoadRange.append(float(a[0]))
    LoadCycles.append(float(a[1]))


data={
      "LoadRange":[30.0,24.0,24.0,22.0,22.0,21.0,20.0,20.0,18.0,18.0,16.0,16.0,16.0,15.0,14.0,14.0,13.0,12.0,12.0,11.0,11.0,10.0,10.0,9.5,9.0,9.0,8.0,8.0,8.0,7.0,7.0,6.0,6.0,6.0,5.8],
      "LoadCycles":[3.462,6.239,7.063,10.189,10.307,11.632,12.593,16.732,20.029,22.292,25.587,27.244,45.298,45.496,66.416,73.128,90.897,95.377,118.807,119.701,192.865,241.237,336.543,345.543,353.876,528.773,886.545,948.345,1021.174,1413.543,1767.506,2081.14,2530.311,2561.642,2692.952]
      }

new_data=pd.DataFrame(data)
 

def fatigue(new_data,which):
    new_data=pd.DataFrame(data)
    
    if which=='plot':
        x=new_data['LoadCycles']
        y=new_data['LoadRange']
        plt.title('F-N',fontsize=22,c='#1f77b4')
        plt.plot(x,y,'o-r',c='#9467bd')
        plt.xlabel('LoadCycles(n)',fontsize=20,c='#2ca02c')
        plt.ylabel('LoadRange(KN)',fontsize=20,c='#bcbd22')
        plt.xlim(0,2700)
        plt.ylim(0,40)
        plt.show()
    
    if which=='max_LoadRange':
        m=new_data['LoadRange'].max()
        return m
    
    if which=='max_LoadCycles':
        n=new_data['LoadCycles'].max()
        return n
    
   
adrs='C:\\Users\\Brooz\\Desktop\\datafatigue.csv'    
fatigue(adrs,'plot')
fatigue(adrs,'max_LoadCycles')
fatigue(adrs,'max_LoadRange')





