#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 20:10:45 2025

@author: apm

functions


"""

#import-----------------------------------------
import math
import statistics
import cmath
import random
import numpy as np
import pandas as pd
import matplotlib as plt






#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
#==================================================================================
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

