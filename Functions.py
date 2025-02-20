'''
Functions.py :









'''

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


#import-----------------------------------------
import math
import statistics
import cmath
import random
import numpy as np
import pandas as pd
import matplotlib as plt





def Activation_Energy(k,k0,T):
    K=math.log(k)
    R=8.3144598

    K0=math.log(k0)
    Ea=(K0-K)*R*T
    return Ea 


def Atomic_Packing_Factor(radius , crystal_structure):
    
    '''
    Parameters
    ----------
    
    crystal_structure: Type of the crystal (Fccc , Bcc or ...)
    radiuse (float): Atomic Atomic radius in crystal structure (nm)
    
    Returns:
    Atomic packing factor for SC , BCC , FCC
    
    '''

    if crystal_structure.upper() == "SC":
        volume = (4/3) * math.pi * radius ** 3
        cell_volume = radius * 2 
        APF = volume / cell_volume ** 3
        
    elif crystal_structure.upper() == "BCC":
        volume = (4/3) * math.pi * radius ** 3
        cell_volume = 4 * radius * math.sqrt(3) / 2
        APF = volume / cell_volume ** 3 
        
    elif crystal_structure.upper() == "FCC":
        volume = (4/3) * math.pi * radius ** 3
        cell_volume = 16 * radius / math.sqrt(2)
        APF = volume / cell_volume ** 3
        
    else:
        print("Invalid crystal structure. Please enter SC , BCC or FCC")
        return None
    
    return APF








def Activity_Coef(wB,wC,wD,eBB,eCB,eDB):
    '''
           
    Parameters
    ----------
    wB  : float
          Weight percent of B     
    wC  : float
          Weight percent of C       
    wD  : float
          Weight percent of D       
    eBB : float
          Interaction coefficient of B and B       
    eCB : float
          Interaction coefficient of C and B       
    eDB : float
          Interaction coefficient of D and B       
    
    Returns
    -------
    fB : float
          Activity coefficient of B

    '''
    
    fB=e**(eBB*wB+eCB*wC+eDB*wD)
    
    return fB

def Arithmetic_Sequence(start_num,common_difference,n):
    '''
An arithmetic sequence is an ordered set of numbers that have a common difference between each consecutive term.
    Parameters
    ----------
    start_num : int
        the first term in the sequence.
    common_difference : int
        the common difference between terms.
    n : int
        number of terms.

    Returns
    -------
    a_n	:int
    	the nᵗʰ term in the sequence

    '''
    a_n = start_num + ((n - 1)*common_difference)
    return a_n

def Aeroscope_Stress_Concentration(max_stress,nominal_stresss):
    K=max_stress/nominal_stresss
    return K







def BMI_Calculation(W,H):
    '''
    
    This function calculates body mass index
    Parameters
    ----------
    W : float
        Weight in kilograms
    H : float
        height in meters

    Returns
    float
    BMI_Calculation

    '''
    return(W/H**2)


def Biomaterial_Degredation_Rate(W1,W2,T):
    '''
    This function calculates the degradation  rate of biomaterials
    
    Parameters
    ----------
    W1 : int
        initial mass
    W2 : int
        final mass
    T : int
        time of degredation

    Returns
    int
   Biomaterial_Degredation_Rate

    '''
    return((W1-W2)/W1*100/T)     


def Copolymer_Type(Copolymer):
    co_list = Copolymer.split()
    change = []
    for i in range(1,len(co_list)):
        if co_list[i] != co_list[i-1]:
            change.append(True)
        else:
            change.append(False)

    if change.count(False) == 0:
        print('alternative')
    elif change.count(True) == 0:
        print('This is not a copolymer')
    elif change.count(True)==1:
        print('Block')
    else: 
        print('Random')
        
def Corrosion_Rate(W,D,A,t):
    '''
    

    Parameters
    ----------
    W : int or float
        The change in weight of specimen (mg)      
    D : int or float
        The density of specimen (g/cm^3)      
    A : int or float
        The surface area of specimen (in^2)      
    t : int or float
        The time(h)

    Returns
    -------
    CR : int or float
         The corrosion rate (mpy)

    '''
    
    CR=534.6*W/(D*A*t)
    return CR


def Circle_Area(radius):
    '''
    

    Parameters
    ----------
    radius : int
        radius of circle.

    Returns
    -------
    circle_area:int
        area of circle.

    '''
    circle_area=(radius**2)*pi
    return circle_area



    
def Circle_Perimeter(radius):
    '''
    

    Parameters
    ----------
    radius : int
        radius of circle.

    Returns
    -------
    circle_perimeter: int
        perimeter of circle.

    '''
    circle_perimeter=2*pi*radius
    return circle_perimeter
    


def Circle_Area(r):
    x=r**2*pi
    return x


def Circle_Surrondings(r):
    x=r*ta
    return x



def Density(m, V):
    '''
    It's the formula for obtaining density.
    
    Parameters:
    ----------
    m : float
        Mass
    
    V: float
        volume
    '''
    den = m / V       
    return den

def Drag_Force(Velocity,Fluid_Coefficent,Fluid_Density,cross_sectional_area,):
   D = ((Velocity**2)*(Fluid_Density)*(Fluid_Coefficent)*(cross_sectional_area))/2
   return D
    

def Density (mass,volume):
    d=mass/volume
    return d


def Error_Function(z):
    '''
    

    Parameters
    ----------
    z : int or float
        Error function argument

    Returns
    -------
    erf : float
        error function value

    '''
    erf=0
    if z<0:
        z=-z
    t=0
    d=0.00001
    while t<z:
        f1=e**(-t**2)
        f2=e**(-(t+d)**2)
        erf=erf+(2/(pi**0.5))*((f1+f2)*d/2)
        t=t+d
        
    erf=int(erf*1000000)/1000000
    return erf
     


def Encapsulation_Efficiency(W1,W2 ):
    '''
This function calculates the percentage of drug loaded in the carrier during drug delivery

    Parameters
    ----------
    W1 : float
        weight of free drug
    W2 : float
        weight of total drug

    Returns
    float
    Encapsulation_Efficiency
        

    '''
                
    return 1-(W1/W2)*100


def First_Row_Pascal_Triangle(k):
    
    result=1
    
    for i in range(1,k+1):
        result*=i
        
    return result


    
def Fibonachi_Sequence (N):
    fibo = [0,1]
    first = 0
    second = 1
    for i in range(2,N):
        new = first + second
        fibo.append(new)
        first = second
        second = new
    return fibo

def Factorial(a):
    '''
    The product of all positive integers less than or equal to a given positive integer.

    Parameters
    ----------
    a : int
       

    Returns
    -------
    factorial: int
      the product of all positive integers less than or equal to a given positive integer .

    '''
    if a==0:
        return 1
    else:
        i=1
        factorial=a
        while a!=i:
            factorial=factorial*(a-i)
            i=i+1
        return factorial


def Fabric_GSM(Warp,Weft,Warp_Count_Nm,Weft_Count_Nm,Shrinkage_Percent=5):
   '''
    This function calculates weight fabric in GSM unit.

    Parameters
    ----------
    Warp : int or float
        The number of warps in 1 cm.
    Weft : int or float
        The number of wefts in 1 cm.
    Warp_Count_Nm : int or float
        Warp yarn metric count.
    Weft_Count_Nm : int or float
        Weft yarn metric count.
    Shrinkage_Percent : int or float, optional
        The percentage difference in the length of woven and non-woven yarn. The default is 5.
    Fabric_GSM : int or float
        Result.

    '''
   Fabric_weight= ((Warp*100)/Warp_Count_Nm )+ ((Weft*100)/Weft_Count_Nm )
   Fabric_GSM= Fabric_weight * (1+(Shrinkage_Percent/100))
   return Fabric_GSM

def Fabric_Drape_Coefficient(fabric_weight,fabric_thickness,bending_length):
    '''
    This function estimates the drape coefficient of fabric according to 3 factors:

    Parameters
    ----------
    fabric_weight : int or float
        weight of fabric.  
    fabric_thickness : int or float
        thickness of fabric. 
    bending_length : int or float
        the length difference. 
    Drape_Coefficient : int or float
        Result.

    '''
    Drape_Coefficient = (fabric_weight*bending_length)/(fabric_thickness**2)
    return Drape_Coefficient


def Fabric_Porosity(air_volume,total_volume):
    '''
    

    Parameters
    ----------
    air_volume : int
        enter the volume of pores in fabric in mm^3.
    total_volume : int
        Enter the total volume of fabric in mm^3.
    Returns
    -------
    Fabric_porosity : int
        This function calculates the fabric porosity in mm^3.

    '''
    FP=total_volume-air_volume
    return FP


def Fabric_weight(density,area):
    '''
    

    Parameters
    ----------
    density : int
        enter the density of fabric in g/mm^2.
    area : int
        Enter the area of fabric in mm^2.
    Returns
    -------
    Fabric_weight : int
        This function calculates the fabric weight in g.

    '''
    FW=density*area
    return FW


"""
# 3
# calculate Fibonachi_Sequence by get N
"""
def Fibonachi_Sequence (N):
    fibo = 0
    for i in range(0,N+1):
        fibo = fibo + i
    return fibo


def Geometric_Sequence(first_variable,second_variable):
    '''
    This function obtains the result of a geometric progression

    Parameters
    ----------
    first_variable : int
        The value of the first variable.
    second_variable : int
       The value of the second variable.

    Returns
    -------
    int
        gives the result of the geometric exponential function.

    '''
    if second_variable>1:
       m=(first_variable*(Geometric_Sequence(first_variable,second_variable-1)))/2 
       
    else:
        
        return 1
    return m 



def Heat_Transfer_Rate (Thermal_conductivity,area, tempreture_difference):
    Heat_Transfer_Rate=Thermal_conductivity*area*tempreture_difference
    return Heat_Transfer_Rate
'This function is used to calculate heat transfer rate by using thermal conductivity, area of transformation, and tempreture of to sides of transformation.'
'In this formula, heat transfer rate is in Btu/(hr*square ft*F), area is in square ft, and tempreture difference is in F.'


  
def Half_Inhibitory_Concentration(A1,A2,A3):
    '''
This function determines the drug’s efficacy in biological process inhibition     

    Parameters
    ----------
    A1 : int
        absorbance of the experimental wells 
    A2 : int
        absorbance of the control wells
    A3 : int
        absorbance of the blank wells

    Returns
     float
        Half _Inhibitory_Concentration 

    '''
    return 1-(A1-A3/A2-A3)*100
  
def Hooke(strain,young_modulus):
    stress=young_modulus*strain
    return stress







def Ideal_Gas_low(R = "(L.atm) / (K.mol)", V = 1, n = 1, T = 0):
    
    '''
    It calculate the Pressure of ideal gas.
    
    Parameters:
    ----------
    P : float
        Pressure
        
    V : float
        The amount of Volume 
        
    n : float 
    
    The amount of substance
    
    T : float
        Temperature
    
    R : float
        ideal gas constant
        
    '''
    # V in Litr, Pressure in atm, n in mol.
    if R == "(L * atm) / (K * mol)":
        P = (n * T * 0.082057) / V
        return P
    
    if R == " J / (K.mol)":
        P = (n * T * 8.314472) / V
        return P
    
    if R == "((m**3).atm) / (K.mol)":
        P = (n * T * 8.205745 * (10 **(-5))) / V
        return P
    
    if R == "(L.kPa) / (K.mol)":
        P = (n * T * 8.314472) / V
        return P
    
    if R == "(((m**3).Pa) / (K.mol))":
        P = (n * T * 8.314472) / V
        return P
    
    if R == "((cm ** 3).atm) / (K.mol)":
        P = (n * T * 82.05745) / V
        return P
    
    if R == "(L.mbar) / (K.mol)":
        P = (n * T * 83.14472) / V
        return P
    
    if R == "((m**3).bar) / (K.mol)":
        P = (n * T * (8.314472 * (10**(-5)))) / V
        return P
    



def Mandaliof_Properties(E):
    '''
    

    Parameters
    ----------
    E : string
        The chemical symbol of the Element

    Returns
    -------
    Nothing

    '''
    PT=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
        'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr',
        'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
        'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']
    loc=PT.index(E)
    Column=[1,18,1,2,13,14,15,16,17,18,1,2,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
        15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,3,3,3,3,3,3,3,
        3,3,3,3,3,3,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,3,3,3,3,3,3,3,
        3,3,3,3,3,3,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
    Row=[1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,
         5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,
         7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7]
    R=[31,28,128,112,84,76,71,66,57,58,166,141,121,111,107,105,102,106,203,176,170,160,153,139,139,132,152,124,132,122,122,120,
        119,120,120,116,220,195,190,175,164,154,147,146,142,139,145,144,142,139,139,138,139,140,244,215,207,204,203,201,199,198,198,196,
        194,192,192,189,190,187,187,175,170,162,151,144,141,136,136,132,145,146,148,140,150,150,260,221,215,206,200,186,190,187,180,169,
        168,'N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A']
    A_W=[1.008,4.003,6.94,9.012,10.81,12.01,14.007,15.999,18.998,20.180,22.990,24.305,26.981,28.085,30.974,32.06,35.45,39.948,39.098,40.078,44.956,47.867,50.941,51.996,54.938,55.845,58.933,58.693,63.546,65.38,69.723,72.63,
        74.922,78.971,79.904,83.798,85.468,87.62,88.906,91.224,92.906,95.95,98,101.07,102.906,106.42,107.868,112.414,114.818,118.71,121.76,127.6,126.904,131.293,132.905,137.327,138.905,140.116,140.908,144.242,145,150.36,151.964,157.25,
        158.925,162.5,164.930,167.259,168.934,173.045,174.967,178.49,180.948,183.84,186.207,190.23,192.217,195.084,196.967,200.592,204.38,207.2,208.980,209,210,222,223,226,227,232.038,231.036,238.029,237,244,243,247,
        247,251,252,257,258,259,262,267,270,269,270,270,278,281,281,285,286,289,289,293,293,294]
    T_m=[-259.14,'N/A',180.54,1287,2075,3550,-210.1,-218.3,-219.6,-248.59,97.72,650,660.32,1414,44.2,115.21,-101.5,-189.3,63.38,842,1541,1668,1910,1907,1246,1538,1495,1455,1084.62,419.53,29.76,938.3,
        817,221,-7.3,-157.36,39.31,777,1526,1855,2477,2623,2157,2334,1961,1554.9,961.78,321.07,156.6,231.93,630.63,449.51,113.7,-111.8,28.44,727,920,798,931,1021,1100,1072,822,1313,
        1356,1412,1474,1497,1545,819,1663,2233,3017,3422,3186,3033,2466,1768.3,1064.18,-38.83,304,327.46,271.3,254,302,-71,20.9,700,1050,1750,1572,1135,644,640,1176,1345,
        1050,900,860,1500,830,830,1600,'N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A']
    T_b=[-252.87,-268.93,1342,2470,4000,4027,-195.79,-182.9,-188.12,-246.08,883,1090,2519,2900,280.5,444.72,-34.04,-185.8,759,1484,2830,3287,3407,2671,2061,2861,2927,2913,2562,907,2204,2820,
        614,685,59,-153.22,688,1382,3345,4409,4744,4639,4265,4150,3695,2963,2162,767,2072,2602,1587,988,184.3,-108,671,1870,3464,3360,3290,3100,3000,1803,1527,3250,
        3230,2567,2700,2868,1950,1196,3402,4603,5458,5555,5596,5012,4428,3825,2856,356.73,1473,1749,1564,962,350,-61.7,650,1737,3200,4820,4000,3900,4000,3230,2011,3110,
        'N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A','N/A']
    Name=['Hydrogen','Helium','Lithium','Beryllium','Boron','Carbon','Nitrogen','Oxygen','Fluorine','Neon','Sodium','Magnesium','Aluminum','Silicon','Phosphorus','Sulfur','Chlorine','Argon','Potassium','Calcium','Scandium','Titanium','Vanadium','Chromium','Manganese','Iron','Cobalt','Nickel','Copper','Zinc','Gallium','Germanium',
        'Arsenic','Selenium','Bromine','Krypton','Rubidium','Strontium','Yttrium','Zirconium','Niobium','Molybdenum','Technetium','Ruthenium','Rhodium','Palladium','Silver','Cadmium','Indium','Tin','Antimony','Tellurium','Iodine','Xenon','Caesium','Barium','Lanthanum','Cerium','Praseodymium','Neodymium','Promethium','Samarium','Europium','Gadolinium',
        'Terbium','Dysprosium','Holmium','Erbium','Thulium','Ytterbium','Lutetium','Hafnium','Tantalum','Tungsten','Rhenium','Osmium','Iridium','Platinum','Gold','Mercury','Thallium','Lead','Bismuth','Polonium','Astatine','Radon','Francium','Radium','Actinium','Thorium','Protactinium','Uranium','Neptunium','Plutonium','Americium','Curium',
        'Berkelium','Californium','Einsteinium','Fermium','Mendelevium','Nobelium','Lawrencium','Rutherfordium','Dubnium','Seaborgium','Bohrium','Hassium','Meitnerium','Darmstadtium','Roentgenium','Copernicium','Nihonium','Flerovium','Moscovium','Livermorium','Tennessine','Oganesson']
    
    if loc>95:
        P='Solid'
    elif loc==1:
        P='Gas'
    else:
        if 25<T_m[loc]:
            P='Solid'
        elif T_m[loc]<25<T_b[loc]:
            P='Liquid'
        else:
            P='Gas' 
    print('The properties of the element you are looking for are:\n\n\nFull name:                  ',Name[loc],'\nLocation in Periodic Table:  Row=',Row[loc],'  Column=',Column[loc],
          '\nAtomic number:              ',loc+1,'\nAtomic weight:              ',A_W[loc],' g/mol\nAtomic radius:              ',R[loc],' pm\nMelting temperature:        ',T_m[loc],end='')
    if loc==32:
        print(' (Above 28 atm)',end='')
    print('\nBoiling temperature:        ',T_b[loc],' \nStandard phase:             ',P)


def Mtt_Test(C1,C2):
    '''
 This function measures the metabolic activity of the cell
    
    Parameters
    ----------
    C1 : int
        viable cells 
    C2 : int
        total cells

    Returns
    float
        Mtt_Test

    '''
        
    
    return (C1/C2)*100
    
    
def Nanoparticle_Surface_Area(Shape,Diameter=0,a=0):
    """
    Calculating the surface area of nanoparticle by determinig the shape
    
    Parameters
    ----------
    Shape:str
    (sphere,cube or tetrahedron)
    
    Diameter : int
    DESCRIPTION. The default is 0.
    (it is needed for spherical nanoparticles)
    unit(nm)
       
    a : int
    dimention
    unit(nm)
    DESCRIPTION. The default is 0.
    (it is needed for cubic or tetragonal nanoparticles)
    unit(nm)
    
    Returns
    -------
    Nanoparticle_Surface_Area : int
    unit(nm^2) 

    """
    Shape=input('enter the shape of nanoparticle')
    if Shape=='sphere':
        Diameter=int(input('enter the shape of nanoparticles diameter'))
        Nanoparticle_Surface_Area=Pi*(Diameter**2)
        print(Nanoparticle_Surface_Area)
        return Nanoparticle_Surface_Area
    elif Shape=='cube':
        a=int(input('enter a '))
        Nanoparticle_Surface_Area=6*(a**2)
        print(Nanoparticle_Surface_Area)
        return Nanoparticle_Surface_Area
    elif Shape=='tetrahedron':
        a=int(input('enter a '))
        Nanoparticle_Surface_Area=(3**0.5)*(a**2)
        print(Nanoparticle_Surface_Area)
        return Nanoparticle_Surface_Area
    else: 
        print('please retry and enter the needed parameters corectly')


def Nanoparticle_Aspect_Ratio(lenght=1,width=1,height=1):
    """
    Calculating the Nanoparticle_Aspect_Ratio 

    Parameters
    ----------
    lenght : int
        lenght of a nanoparticle
        DESCRIPTION. The default is 1.
        unit(nm)
    width : int
        width of a nanoparticle
        DESCRIPTION. The default is 1.
        unit(nm)
    height : int
       height of a nanoparticle
       DESCRIPTION. The default is 1.
       unit(nm)

    Returns
    -------
    Nanoparticle_Aspect_Ratio : float
        the ratio of maximum dimention of a nanoparticle to minimum dimention of it.

    """
    lenght=int(input('please enter lenght'))
    width=int(input('please enter width'))
    height=int(input('please enter height'))
    maximum=max(lenght,width,height)
    minimum=min(lenght,width,height)
    Nanoparticle_Aspect_Ratio=maximum/minimum
    print('Nanoparticle_Aspect_Ratio=',Nanoparticle_Aspect_Ratio)
    return Nanoparticle_Aspect_Ratio

def Nanoparticle_Volume(Shape,Diameter=0,a=0):
    """
    Calculating the Volume of a nanoparticle by determinig the shape
    
    Parameters
    ----------
    Shape:str
    (sphere,cube or tetrahedron)
    
    Diameter : int
    DESCRIPTION. The default is 0.
    (it is needed for spherical nanoparticles)
    unit(nm)
       
    a : int
    dimention
    unit(nm)
    DESCRIPTION. The default is 0.
    (it is needed for cubic or tetragonal nanoparticles)
    unit(nm)
    
    Returns
    -------
    Nanoparticle_Volume : int
    unit(nm^3) 

    """
    Shape=input('enter the shape of nanoparticle')
    if Shape=='sphere':
        Diameter=int(input('enter the shape of nanoparticles diameter'))
        Nanoparticle_Volume=Pi*((Diameter**3)/6)
        print(Nanoparticle_Volume)
        return Nanoparticle_Volume
    elif Shape=='cube':
        a=int(input('enter a '))
        Nanoparticle_Volume=a**3
        print(Nanoparticle_Volume)
        return Nanoparticle_Volume
    elif Shape=='tetrahedron':
        a=int(input('enter a '))
        Nanoparticle_Volume=(2**0.5)*(a**3)/12
        print(Nanoparticle_Volume)
        return Nanoparticle_Volume
    else: 
        print('please retry and enter the needed parameters corectly')



def poisson(transverse_strain,axial_strain):
    v=-(transverse_strain/axial_strain)
    return v

def pH_Calculator(a, c, K1=0, K2=0, K3=0):

    
    
    '''
    Parameters
    ----------
    a : int
        If you have a strong acid, enter 1. For strong base, enter 2. For weak acid, enter 3. For weak base, enter 4.
    c : float
        Enter the concentration of acid or base solution (mol/Lit).
    K1 : float
        Enter the first dissociation constant of the weak acid or base. 
        The dissociation constant for weak acids and bases must be lower than one(K<1).
        For polyprotic acid or bases (the acid or bases that have more than one hydrogen or hydroxide ion per molecule), only enter only dissociation coefficients smaller than 1. (K<1)
    K2 : float
        Enter the second dissociation constant of the weak acid or base.
    K3 : float
         Enter the third dissociation constant of the weak acid or base.    
    Returns
    -------
    pH : float
        The pH of acid/base solution.
    '''
        
    
    
    
    
    
    import math as m
    pHw=-m.log10(Kw)/2
    Hw=10**(-pHw)
    H=[0,0,0]
    x=[0,0]
    pH=0
    c= float(c)
    k=[K1,K2,K3] 
    if a == 1:
            H[0]=c
            x[0]=(-(H[0]+k[0])+m.sqrt((H[0]+k[0])**2+4*k[0]))/2
            x[1]=(-(H[0]+k[0])-m.sqrt((H[0]+k[0])**2+4*k[0]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[1]=H[0]+x[i]
            x[0]=(-(H[1]+k[1])+m.sqrt((H[1]+k[1])**2+4*k[1]))/2
            x[1]=(-(H[1]+k[1])-m.sqrt((H[1]+k[1])**2+4*k[1]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[2]=H[1]+x[i]
                    pH = -m.log10(H[2]+Hw)
                    return pH
                 
            
            
    elif a == 2:
            H[0]=c
            x[0]=(-(H[0]+k[0])+m.sqrt((H[0]+k[0])**2+4*k[0]))/2
            x[1]=(-(H[0]+k[0])-m.sqrt((H[0]+k[0])**2+4*k[0]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[1]=H[0]+x[i]
            x[0]=(-(H[1]+k[1])+m.sqrt((H[1]+k[1])**2+4*k[1]))/2
            x[1]=(-(H[1]+k[1])-m.sqrt((H[1]+k[1])**2+4*k[1]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[2]=H[1]+x[i]
                    pH = 14+m.log10(H[2]+Hw)
                    return pH
            
            


    elif a == 3:
        if k[0]==0:
            print('Enter the dissociation constants of the weak acid')
        else:
            H[0]=(-k[0]+m.sqrt(k[0]**2+4*c*k[0]))/2
            x[0]=(-(H[0]+k[1])+m.sqrt((H[0]+k[1])**2+4*k[1]))/2
            x[1]=(-(H[0]+k[1])-m.sqrt((H[0]+k[1])**2+4*k[1]))/2
            for i in range(2):
                if  type(x[i])!='complex' and x[i]>=0:
                    H[1]=H[0]+x[i]
            x[0]=(-(H[1]+k[2])+m.sqrt((H[1]+k[2])**2+4*k[2]))/2
            x[1]=(-(H[1]+k[2])-m.sqrt((H[1]+k[2])**2+4*k[2]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[2]=H[1]+x[i]
                    pH = -m.log10(H[2]+Hw)
                    return pH
              
            


    elif a == 4:
        if k[0]==0:
            print('Enter the dissociation constants of the weak base')
        else:
            H[0]=(-k[0]+m.sqrt(k[0]**2+4*c*k[0]))/2
            x[0]=(-(H[0]+k[1])+m.sqrt((H[0]+k[1])**2+4*k[1]))/2
            x[1]=(-(H[0]+k[1])-m.sqrt((H[0]+k[1])**2+4*k[1]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[1]=H[0]+x[i]
            x[0]=(-(H[1]+k[2])+m.sqrt((H[1]+k[2])**2+4*k[2]))/2
            x[1]=(-(H[1]+k[2])-m.sqrt((H[1]+k[2])**2+4*k[2]))/2
            for i in range(2):
                 if  type(x[i])!='complex' and x[i]>=0:
                    H[2]=H[1]+x[i]
                    pH = 14+m.log10(H[2]+Hw)
                    return pH
              
 






            
  
    
 
    
def pH_Which (pH):
    if pH<7:
        print ('This solution is acidic')
    elif pH ==7:
        print ('This solution is neutral')
    elif pH>7:
        print('This solution is basic')
           
        
        






def Principal_Stress(Sx,Sy,Sz,Txy,Tyz,Txz):
    '''
    make sure stress value is less than 1000       
    Parameters
    ----------
    Sx  : (int or float)
          Normal stress along x       
    Sy  : (int or float)
          Normal stress along y       
    Sz  : (int or float)
          Normal stress along z       
    Txy : (int or float)
          Shear stress along xy       
    Tyz : (int or float)
          Shear stress along yz       
    Txz : (int or float)
          Shear stress along xz       
    Returns
    -------
    S_P : [S1,S2,S3]
          Principal stresses

    '''
    
      
    a=1
    b=Sx+Sy+Sz
    c=Sx*Sy+Sy*Sz+Sx*Sz-Txy**2-Tyz**2-Txz**2
    d=Sx*Sy*Sz+2*Txy*Tyz*Txz-Sx*Tyz**2-Sy*Txz**2-Sz*Txy**2
    
    
    # aS^3-bS^2+cS-d=0

    S_P=[0,0,0]
    #-------------Numerical Calculation---------------
    sp=0
    for i in range(2001):
        x0=-1000+i
        f0=a*x0**3-b*x0**2+c*x0-d
        x1=x0+1
        f1=a*x1**3-b*x1**2+c*x1-d
        if f0>-10 and f0<10:
            S_P[sp]=x0
            if sp==2:
                break
            sp=sp+1
        else:
            if f0*f1<0:
                while i>-1:
                    x=(x0+x1)/2
                    f=a*x**3-b*x**2+c*x-d
                    if f>-10 and f<10:
                        S_P[sp]=x
                        if sp==2:
                            break
                        sp=sp+1
                        break
                    elif f0*f<0:
                        x1=x
                        f1=f
                    else:
                        x0=x
                        f0=f
            else:
                continue
                        
    
    return S_P


def Pythagorean(side1,side2):
    '''
    It should be an orthogonal triangle and this function gives you hypotenuse

    Parameters
    ----------
    side1 : int
        side of "orthogonal triangle".
    side2 : int
        side of "orthogonal triangle".

    Returns
    -------
    hypotenuse: int
       hypotenuse of "orthogonal triangle".

    ''' 
    hypotenuse =((side1**2)+(side2**2))**(1/2)
    return hypotenuse




def Quadratic_Equation(a,b,c):
    '''
    This function find "x" in equation ax^2 + bx + c = 0

    Parameters
    ----------
    a : int
        known number.
    b : int
        known number.
    c : int
        known number.

    Returns
    -------
   The Quadratic Equation mostly have two answers!
    x1 ,x2 : int
       the unknown.

    '''
    if a==0:
        print('Error: This is not a quadratic equation')
    elif ((b**2)-(4*a*c))<0:
        print("Error:this equation doesn't have an answer")
    else:
        x1=(-b+(((b**2)-(4*a*c))**(1/2)))/(2*a)
        x2=(-b-(((b**2)-(4*a*c))**(1/2)))/(2*a)
        if x1==x2:
            return x1
        else:
            return x1,x2
        
def Rectangle_Area(length,width):
    '''
    This function calculate the area of square too!

    Parameters
    ----------
    length : int
        length of rectangle.
    width : int
        width of rectangle.

    Returns
    -------
    rectangle_area: int
      area of rectangle.

    '''
    rectangle_area=length*width
    return rectangle_area
    
    


def Rectangle_Perimeter(length,width):
    '''
    This function calculate the perimeter of square too!  

    Parameters
    ----------
    length : int
        length of rectangle.
    width : int
        width of rectangle.

    Returns
    -------
    rectangle_perimeter: int
       

    '''
    rectangle_perimeter=2*(length+width)
    return rectangle_perimeter


def root_degree2(a,b,c):
    delta=b**2-(4*a*c)

    if delta>0:
        x=delta*0.5+b/2*a and delta*0.5-b/2*a 
        return x
    if delta==0:
        x=-b/2*a 
        return x
    if delta<0:
        print("no answer")
        
        
def shearrate(Q,p,r,n):
    y=(Q/p*r**3)*(3+1/n)
    return y
#2-formula of shear rate
#F=applied force,A=cross-sectional area,T=shear stress
def shearstress (F,A):
    T=F/A
    return T


def Surface_Area_To_Volume_Ratio(Shape,Diameter=0,a=0):
    """
    Calculating the ratio of nanoparticle's surface to its volume by determinig the shape
            
    Parameters
    ----------
    Shape:str
    (sphere,cube or tetrahedron)
            
    Diameter : int
    DESCRIPTION. The default is 0.
    (it is needed for spherical nanoparticles)
    unit(nm)
               
    a : int
    dimention
    unit(nm)
    DESCRIPTION. The default is 0.
    (it is needed for cubic or tetragonal nanoparticles)
    unit(nm)
            
    Returns
    -------
    Surface_Area_To_Volumr_Ratio : float
    unit(1/nm) 

    """
    Shape=input('enter the shape of nanoparticle')
    if Shape=='sphere':
        Diameter=int(input('enter the shape of nanoparticles diameter'))
        Nanoparticle_Surface_Area=Pi*(Diameter**2)
        Nanoparticle_Volume=Pi*((Diameter**3)/6)
    elif Shape=='cube':
        a=int(input('enter a '))
        Nanoparticle_Surface_Area=6*(a**2)
        Nanoparticle_Volume=a**3
    elif Shape=='tetrahedron':
        a=int(input('enter a '))
        Nanoparticle_Surface_Area=(3**0.5)*(a**2)
        Nanoparticle_Volume=(2**0.5)*(a**3)/12
    else: 
        print('please retry and enter the needed parameters corectly')
    Surface_Area_To_Volume_Ratio=Nanoparticle_Surface_Area/Nanoparticle_Volume
    print( Surface_Area_To_Volume_Ratio)
    return  Surface_Area_To_Volume_Ratio




def Sphere_Volume (R):
    V = (4/3)* math.pi * R**3
    return V

#print("Sphere_Volume = ",Sphere_Volume(2))

"""
# 2
# calculate Area of Sphere by get R
"""
def Sphere_Area (R):
    S = 4 * math.pi * R**2
    return S

#print("Sphere_Area = ",Sphere_Area(2))


def Specific_Impulse(Thrust,Propellent_flowrate):
    d=Thrust/(Propellent_flowrate*g0)
    return d


def Sphere_Volume(radius):
    '''
    

    Parameters
    ----------
    radius : int
        radius of sphere.

    Returns
    -------
    sphere_volume: int
        volume of sphere.

    '''
    sphere_volume=(4*(radius**3)*pi)/3
    return sphere_volume





    
    
def Sphere_Surface_Area(radius):
    '''
    

    Parameters
    ----------
    radius : int
        radius of sphere.

    Returns
    -------
    sphere_surface_area: int
        surface area of sphere .

    '''
    sphere_surface_area=(4*pi*(radius**2))
    return sphere_surface_area
    
    
    
def Sample_Size_Calculation(N,e):
    '''
    
This  function estimates the number of evaluable subjects required for achieving desired statistical significance for a given hypothesis. 

    Parameters
    ----------
    N : int
        population of the study 
    e : float
        margin error

    Returns
    int
      Sample_Size_Calculation

    '''
    return((N)/(1+N)*(e**2))
    
    
      
#PART.2.8



def Stress(strain,young_modulus):
    '''
    #The Young's modulus (E) is a property of the material that tells us how easily it can stretch and deform 
    and is defined as the ratio of tensile stress (σ=S) to tensile strain (ε=Ts)

    Parameters
    ----------
    strain : tensile strain (ε)
    
    young_modulus : Modulus of Elasticity (N/m2)

    Returns
    -------
    S : tensile stress (σ)

    '''
       
    global S
    S=strain*young_modulus
    return S






def Trapezium_Area(Base1,Base2,Height):

    '''
	Parameters
    ----------
    Base1 : int
       base of trapezium
    
    Base2 : int
       base of trapezium
      
    Height : int
       height of trapezium

    Returns
    -------
    Area : int
        Area of trapezium
    '''
    Area=((Base1+Base2)*Height)/2
    return Area





def Trapezium_Perimeter(Base1,Base2,Side1,Side2):
    '''
    Parameters
    ----------
    Base1 : int
       base of trapezium.
    
    Base2 : int
        base of trapezium.
    
    Side1 : int
        side of trapezium.
    
    Side2 : int
        side of trapezium.

    Returns
    -------
    Perimeter : int
       perimeter of trapezium

    '''
    Perimeter=Base1+Base2+Side1+Side2
    return Perimeter



def Tensile_Strength (f,a):
    
    '''
    Parameters
    ----------
    f : int
       force required to break
       
    a : int
       cross sectional area
    -------
    c : int
       Tensile_Strength
    '''
    c=f/a
    return c


#tensile strength
def tensilestrength(Maxloud,cross_sectional_area):
    tensilestrength=Maxloud/cross_sectional_area
    return tensilestrength
    


def Triangle_Environment (first_side,second_side,third_side):
    '''
    This function gets the perimeter of the triangle

    Parameters
    ----------
    first_side : float
        The size of the first side.
    second_side : float
       The size of the second side.
    third_side : float
        The size of the third side.

    Returns
    -------
    Environment : float
       The output is the perimeter of the triangle.

    '''
    
    Environment=first_side+second_side+third_side
    return Environment
   


def Triangle_Area(Height,rule):
    '''
    This function obtains the area of ​​the triangle

    Parameters
    ----------
    Height : float
        The size of the height of the triangle.
    rule : float
        The size of the base of the triangle is.

    Returns
    -------
    area : float
        The output is the area of ​​the triangle.

    '''
    area=(Height*rule)/2
    return (area) 


#دقیقا مطمِن نبودم که ورودی اول چیه و 
#با توجه به فرمول هایی که پیدا کردم احتمال دادم مقاومت حرارتی باشه 
#اما اگر اشتباه متوجه شدم ممنون میشوم بهم اطلاع دهید تا درستش کنم.

def Thermal_Conductivity (thermal_conductance,thickness) :
    R=thermal_conductance
    L=thickness
    return L/R
'''
    تابع 
thermal conductivity یا K   ضریب هدایت حرارتی یا رسانندگی گرمایی
   وات بر کلوین متر (W·K−۱·m−۱)
را به ما میدهد 

متغیر های تابع به صورت زیر اند
thermal_conductance یا مقاومت حرارتی یا R 
وات بر کلوین متر2 ((W·K−1·m−2))

thickness یا ضخامت یا L 
متر (m)
'''



#3-Formula of viscosity
#n=viscosity
def viscosity (y,T):
    n=T/y
    return n



def Young_Modulus(stress,strain):
    E=stress/strain
    return E





#---CHATGPT
def ohms_law(voltage, current, resistance):
    return voltage - current * resistance


def newtons_second_law(force, mass, acceleration):
    return force - mass * acceleration




def kinetic_energy(mass, velocity):
    return 0.5 * mass * velocity**2


def coulombs_law(charge1, charge2, distance):
    k_constant = 8.9875e9  # Coulomb's constant
    return k_constant * (charge1 * charge2) / distance**2

def pythagorean_theorem(a, b):
    return (a**2 + b**2)**0.5

def boyles_law(initial_volume, initial_pressure, final_volume):
    return initial_pressure * initial_volume - final_volume


def ideal_gas_law(pressure, volume, temperature):
    gas_constant = 8.314  # Ideal gas constant
    return pressure * volume - gas_constant * temperature

def youngs_modulus(stress, strain):
    return stress - strain

def heat_transfer(thermal_conductivity, area, temperature_difference, thickness):
    return thermal_conductivity * area * temperature_difference / thickness


def mass_energy_equivalence(mass):
    speed_of_light = 3e8  # Speed of light in m/s
    return mass * speed_of_light**2


def hookes_law(spring_constant, displacement):
    return spring_constant * displacement


def faradays_law(induced_emf, time, magnetic_flux):
    return induced_emf - time * magnetic_flux

def wavelength_frequency_relation(speed_of_light, wavelength, frequency):
    return speed_of_light - wavelength * frequency


def archimedes_principle(density_fluid, volume_displaced, gravitational_acceleration):
    return density_fluid * volume_displaced * gravitational_acceleration

def ideal_diode_equation(current, saturation_current, thermal_voltage):
    return current - saturation_current * (math.exp(current / thermal_voltage) - 1)


def gauss_law(electric_field, surface_area, electric_flux):
    return electric_field * surface_area - electric_flux

def darcys_law(flow_rate, permeability, area, pressure_difference):
    return permeability * area * pressure_difference - flow_rate

def schrodinger_equation(energy, potential_energy, kinetic_energy):
    return energy - (potential_energy + kinetic_energy)

def lorentz_force(charge, velocity, magnetic_field):
    return charge * (velocity.cross(magnetic_field))

def boltzmann_distribution(energy, temperature, boltzmann_constant):
    return math.exp(-energy / (boltzmann_constant * temperature))

def maxwells_equations(electric_field, magnetic_field, charge_density, current_density):
    return electric_field.div() - charge_density, magnetic_field.curl() - current_density

def photoelectric_effect(kinetic_energy, photon_energy, work_function):
    return kinetic_energy - (photon_energy - work_function)

def entropy_change(heat_transfer, temperature):
    return heat_transfer / temperature

def fourier_transform(signal):
    return scipy.fft.fft(signal)


def rayleigh_scattering(intensity, wavelength, particle_size):
    return intensity - (particle_size / wavelength)**4

def nernst_equation(ion_concentration, temperature, faraday_constant, standard_potential):
    return standard_potential - (R * temperature / (n * faraday_constant)) * math.log10(ion_concentration)

def hadamard_product(matrix1, matrix2):
    return np.multiply(matrix1, matrix2)

def elastic_potential_energy(spring_constant, displacement):
    return 0.5 * spring_constant * displacement**2

def rydberg_formula(wavelength, rydberg_constant, principal_quantum_number):
    return 1 / wavelength - rydberg_constant * principal_quantum_number**2

def uncertainty_principle(delta_position, delta_momentum, hbar):
    return delta_position * delta_momentum - hbar / 2




def doppler_effect(observed_frequency, source_frequency, velocity_observer, velocity_source, speed_of_sound):
    return observed_frequency - source_frequency * ((speed_of_sound + velocity_observer) / (speed_of_sound - velocity_source))

def curies_law(magnetic_susceptibility, magnetic_field, temperature):
    return magnetic_susceptibility * magnetic_field - (C / temperature)



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
    

    
class Hydrogen_Energy_Level_Calculater:   #defined the class 
    def Calculate_Energy_Level(sef , n):  #This method is defined with input in n to calculate the energy level
        energy = -13.6 / (n ** 2)         #The energy level formula
        return energy                     #Energy value is returned as the output of the function
    
    def validate_input(self , n):         #This method gets input n  wich represent th hydrogen energy number and check if it's valide or not
        if not isinstance(n , int) or n <= 0:
            raise ValueError("Invalid value!")
            

    
    



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



ef Diffusion_Coefficient(T):
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


def Incompressible_Fluids_Pressure(d,h):  
    P_0=101325#(Pa)
    g=9.81#(m/s^2)
    P=d*g*h+P_0
    return P





def Bouyancy_Force(d,V):
    g=9.81#(m/s^2)
    F_b=d*g*V
    return F_b



def Reynolds_Number_Pipe(d,u,vis,D_H):
    Re=d*u*D_H/vis
    return Re

#functions
def Critical_radius(Temp_difference):
   return (2 * Gama_sl * T_m)/(Latent_heat * Temp_difference)

def Activation_energy(Temp_difference):
    return (16 * Pi * (Gama_sl ** 3) * (T_m ** 2))/(3 * (Latent_heat ** 2) * (Temp_difference ** 2))

def Number_of_unit_cells(Radius):
    return (4 * Pi * (Radius ** 3))/(3 * (converted_A ** 3))

def Nucleation(Temprature, activation):
    return (N_t * math.exp(-(activation/(K * Temprature))))



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



