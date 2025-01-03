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






#-------convertors-----

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


def Area_Converter1(Square_Metre,):
   Square_Cm = Square_Metre*10000
   'this function converts Square_Metre to Square_Cm '
   return Square_Cm



def Area_Converter2(Square_Cm):
   Square_Meter = Square_Cm/10000
   'this function converts Square_Cm to Square_Metre '
   return Square_Meter

def Convertor1 (G_per_Cm3):
    Kg_per_Meter3=G_per_Cm3 * 1000
    return Kg_per_Meter3
 
    
def Convertor2 (Kg_per_Meter3):
    G_per_Cm3=Kg_per_Meter3 / 1000
    return G_per_Cm3

'''
'این دو تابع به عنوان یک convertor  عمل میکنند
اولی برای تبدیل واحد
( گرم بر سانتی متر مکعب 
به 
کیلوگرم بر مترمکعب)
است و دومی برعکس ان را انجام میدهد.

'''
#baraxesh nis
def CelsiusToKelvin(t):
    T=t+273
    return T



def centimeter_to_meter(size):
    '''
    This function converts centimeters to meters

    Parameters
    ----------
    size : int
        Enter an int number to convert centimeters to meters.

    Returns
    -------
    m : int
        The output is the number converted to meters.

    '''
    m=size/100
    return m
    
def mete_to_centimeter (size):
    '''
    This function converts meters  to centimeters

    Parameters
    ----------
    size : int
       Enter an int number to convert meters to centimeters .

    Returns
    -------
    c : int
        The output is the number converted to centimeters.

    '''
    c=size*100
    return c

def Square_Meter_To_Square_Cm(b):
    
    '''
    Parameters
    ----------
    b: int
        Square_meter 
    -------
    c : int
         Square_Cm
    '''
    c =b*10000
    return c

def Square_Cm_To_Square_meter(a):
    
    '''
    Parameters
    ----------
    a : int
        Square_Cm
    -------
    c : int
       Square_Meter
    '''
    c=a/10000
    return c

def Cconverter(c):
    f=(((c*1.8)+32)/32)
    return f


def converter_m_to_mm(meter):
    '''
    

    Parameters
    ----------
    meter : int
        enter the length in meter.
    
    Returns
    -------
    milimeter : int
        This function converts meter into milimeter.

    '''
    milimeter=meter*1000
    return milimeter

def converter_mm_to_m (milimeter):
    '''
    

    Parameters
    ----------
    milimeter : int
        enter the length in milimeter.
    
    Returns
    -------
    meter : int
        This function converts milimeter into meter.

    '''
    meter=milimeter/1000
    return meter







def Cubic_Meter_To_Liter(number_in_Cubic_Meter):
    '''
    This function converts cubic meters to liters.

    Parameters
    ----------
    number_in_Cubic_Meter : int or float
        Number per cubic meter unit. 
    Liter : int or float
        Number per liter unit.

    '''
    Liter= number_in_Cubic_Meter*1000
    return Liter
    

def Liter_To_Cubic_Meter(number_in_Liter):
    '''
    This function converts liters to cubic meters.

    Parameters
    ----------
    number_in_Liter : int or float
        Number per liter unit.
    Cubic_Meter : int or float
        Number per cubic meter unit.

    '''
    Cubic_Meter= number_in_Liter/1000
    return (Cubic_Meter)





# Celcius_To_Kelvin
def Convert_Celcius_To_Kelvin (Celcius):
    Kelvin = Celcius + 273.15
    return Kelvin
'This function is used to convert celcius to kelvin'
'The tempreture in celcius is different from the tempreture in kelvin by 273.15'


#Kelvin_to_celcius
def Convert_Kelvin_to_Celcius (Kelvin):
    Celcius = Kelvin - 273.15
    return Celcius
'This function is used to convert kelvin to celcius'
'The tempreture in celcius is different from the tempreture in kelvin by 273.15'




def Foot_Pound_To_Newton(Foot_Pounds):
    '''
    # This Conventor convert ft-lbs to Nm


    Parameters
    ----------
    Foot_Pound : a unit of torque equal to the force of 1 lb acting perpendicularly to 
    an axis of rotation at a distance of 1 foot.(ft-lbs)

    Returns
    -------
    Newton_Meters : The newton-metre is the unit of torque.(Nm)


    '''
    
      
    global Newton_Meters
    Newton_Meters=Foot_Pounds*1.3558
    return Newton_Meters

def Newton_To_Foot_Pound(Newton_Meters):
    '''
    # This Conventor convert Nm to ft-lbs

    Parameters
    ----------
    Newton_Meters : The newton-metre is the unit of torque.(Nm)

    Returns
    -------
    Foot_Pound : a unit of torque equal to the force of 1 lb acting perpendicularly to 
    an axis of rotation at a distance of 1 foot.(ft-lbs)

    '''    
    
    global Foot_Pound
    Foot_Pound=Newton_Meters*0.7376
    return Foot_Pound


def Fabric_GSM_to_GLM(Fabric_Weight,Fabric_Width):
   '''
    This function converts fabric weight in GSM unit to GLM unit.

     Parameters
     ----------
     Fabric_Weight : int or float
         fabric weight per GSM.
     Fabric_Width : int or float
         width of fabric per inches.
     Fabric_GLM : int or float
        Result.
 
    '''
   Fabric_GLM=(Fabric_Weight*Fabric_Width)/39.37
   return Fabric_GLM
       
def Fconverter(f):
    c=((32*f-32)*(5/9))
    return c



def gpatompa(n):
    Mpa=n*1000
    return Mpa

def mpatogpa(n):
    Gpa=n/1000
    return Gpa




def Gram_To_Picogram(Gram=1):
    """
    converting Gram to Picogram

    Parameters
    ----------
    Gram : float,mass
        DESCRIPTION. The default is 1.

    Returns
    -------
    Picogram : float,mass
       

    """
   
    Gram=int(input ('how many Gram?'))
    Picogram=Gram*1000000000000
    print(Gram,'Gram=',Picogram,'Picogram.')
    return Picogram



def Picogram_To_Gram(Picogram=1):
    """
    converting Picogram to Gram

    Parameters
    ----------
    Picogram : float,mass
        DESCRIPTION. The default is 1.

    Returns
    -------
    Gram : float,mass

    """
    Picogram=int(input ('how many Gram?'))
    Gram=Picogram/1000000000000
    print(Picogram,'Picogram=',Gram,'Gram.')
    return Gram


def Kilogeram_Per_Cubic_Meter_To_Pounds_Per_Cubic_Inch(KgPerCubicMeter):
    L=KgPerCubicMeter*0.0000361273
    return L
def Pounds_Per_Cubic_Inch_To_Kilogeram_Per_Cubic_Meter(LbPerCubicInch):
    Kg=LbPerCubicInch*27679.9
    return Kg




"""
https://abzarek.ir/service-p/length-converter/lightyear-to-Kilometer
# 5
تبدیل کیلومتر به سال نوری
# convert kilometer to Light Year
"""
def KiloMeter_LightYear (km):
    ly = km / 9460730472801.1
    return ly


"""
https://abzarek.ir/service-p/length-converter/lightyear-to-Kilometer
# 4
تبدیل سال نوری به کیلومتر
# convert Light Year to Kilometer 
"""
def LightYear_KiloMeter (ly):
    km = ly * 9460730472801.1
    return km

#print("LightYear_KiloMeter = ",LightYear_KiloMeter(5))




def Micrometer_To_Nanometer(micrometer=1):
    """
    converting micrometer to nanometer 

    Parameters
    ----------
    micrometer : float,dimension
        DESCRIPTION. The default is 1.

    Returns
    -------
    Nanometer : float,dimension
        unit(nm)

    """
    Micrometer=float(input ('how many Micrometer?'))
    Nanometer=Micrometer*1000
    print(Micrometer,'Micrometer=',Nanometer,'Nanometer.')
    return Nanometer


def Nanometer_To_Micrometer(nanometer=1):
    """
    converting nanometer to micrometer

    Parameters
    ----------
    nanometer : float,dimension
      unit (nm)
      DESCRIPTION. The default is 1.
      
    Returns
    -------
    Micrometer : float,dimension
      

    """
    Nanometer=float(input ('how many nanometer?'))
    Micrometer=Nanometer/1000
    print(Nanometer,'nanometer=',Micrometer,'micrometer.')
    return Micrometer

#PART.2.6

def Minute_To_Second (Minute): 
    '''
    This function converts minutes to seconds 

    Parameters
    ----------
    Minute : int
       units of time in minute

    Returns
    
    int
        Minute_to_Second

    '''
       
          
    return (Minute*60)   


def Second_To_Minute (Second):
    '''
This function converts seconds to minutes
        Parameters
    ----------
    Second : int
        units of time in seconds

    Returns
    int
       
      Second_to_Minute
    '''
    
    
    return (Second/60)
        


def Megapascal_To_Pascal(Megapascal):
    '''
    #This Conventor Convert Megapascal to Pascal

    Parameters
    ----------
    Megapascal : 1 Megapascal equals 1,000,000 Pascals.
    

    Returns
    -------
    Pascal : the unit of pressure or stress in SI.
    '''
    
    global Pascal
    Pascal=Megapascal/1000000
    return Pascal

def Pascal_To_Megapascal(Pascal):
    '''
    # This Conventor Convert Pascal to Megapascal

    Parameters
    ----------
    Pascal : the unit of pressure or stress in SI.
    
    Returns
    -------
    Megapascal : 1 Megapascal equals 1,000,000 Pascals.

    '''
    
    global Megapascal
    Megapascal=1000000*Pascal
    return Megapascal





def Newton_TO_Pound_Force(Newton):
     # 1 Pound_Force = 4.448221619 New
     
     
     Pound_Force = Newton / 4.448221619
     '''
     #It converts the Force from Newton to Pound_Force.
     
     Parameters:
     ----------
         
     Newton : float
         Unit musst be newton(N).
         
     '''
     return Pound_Force
 

#_______________


#Pound_Force_TO_Newton    
def Pound_Force_To_Newton(Pound_Force):
    
    Newton = Pound_Force * 4.448221619
    '''
    It converts the Force from Pound_Force to Newton.
    
    Parameters:
    ----------
    
    Pound_Force : float
        Unit musst be lb.
        
    '''
    
    return Newton






def VC_mpsTOkph(mps):
    kph = mps*3.6
    return kph



def VC_kphTOmps(kph):
    mps = kph/3.6
    return mps


def Yarn_Count_Converter(Yarn_Count, Current_System='tex', Desired_System='den'):
    '''
    This function converts yarn count values in different systems.

    Parameters
    ----------
    Yarn_Count : int or float
        Number of yarn count.
    Current_System : str, optional
        Current yarn count system. The default is 'tex'.
    Desired_System : str, optional
        Expected yarn count system. The default is 'den'.
    Yarn_Count : int or float
        Result.

    '''
    sys1=str(Current_System).lower()
    sys2=str(Desired_System).lower()

    if sys1=='tex' and sys2=='dtex':
        Yarn_Count=Yarn_Count*10
        return Yarn_Count
    
    elif sys1=='tex' and sys2=='den':
        Yarn_Count=Yarn_Count*9
        return Yarn_Count

    elif sys1=='tex' and sys2=='nm':
        Yarn_Count=1000/Yarn_Count
        return Yarn_Count
      
    elif sys1=='tex' and sys2=='ne':
        Yarn_Count=590.5/Yarn_Count
        return Yarn_Count
    
    elif sys1=='tex' and sys2=='nw':
        Yarn_Count=885.8/Yarn_Count
        return Yarn_Count
    
    elif sys1=='dtex' and sys2=='tex':
        Yarn_Count=Yarn_Count*0.1
        return Yarn_Count
    
    elif sys1=='dtex' and sys2=='den':
        Yarn_Count=Yarn_Count*0.9
        return Yarn_Count
    
    elif sys1=='dtex' and sys2=='ne':
        Yarn_Count=5905.4/Yarn_Count
        return Yarn_Count
    
    elif sys1=='dtex' and sys2=='nw':
        Yarn_Count=8858/Yarn_Count
        return Yarn_Count
    
    elif sys1=='dtex' and sys2=='nm':
        Yarn_Count=10000/Yarn_Count
        return Yarn_Count
    
    elif sys1=='den' and sys2=='tex':
        Yarn_Count=Yarn_Count/9
        return Yarn_Count
        
    elif sys1=='den' and sys2=='dtex':
        Yarn_Count=Yarn_Count/0.9
        return Yarn_Count
    
    elif sys1=='den' and sys2=='nm':
        Yarn_Count=9000/Yarn_Count
        return Yarn_Count
    
    elif sys1=='den' and sys2=='ne':
        Yarn_Count=5314.9/Yarn_Count
        return Yarn_Count
        
    elif sys1=='den' and sys2=='nw':
        Yarn_Count=7972/Yarn_Count
        return Yarn_Count
    
    elif sys1=='ne' and sys2=='tex':
        Yarn_Count=590.6/Yarn_Count
        return Yarn_Count
    
    elif sys1=='ne' and sys2=='dtex':
        Yarn_Count=5906/Yarn_Count
        return Yarn_Count
    
    elif sys1=='ne' and sys2=='den':
        Yarn_Count=5315/Yarn_Count
        return Yarn_Count
    
    elif sys1=='ne' and sys2=='nm':
        Yarn_Count=1.693*Yarn_Count
        return Yarn_Count
    
    elif sys1=='ne' and sys2=='nw':
        Yarn_Count=1.5*Yarn_Count
        return Yarn_Count
    
    elif sys1=='nm' and sys2=='tex':
        Yarn_Count=1000/Yarn_Count
        return Yarn_Count
    
    elif sys1=='nm' and sys2=='dtex':
        Yarn_Count=10000/Yarn_Count
        return Yarn_Count
    
    elif sys1=='nm' and sys2=='den':
        Yarn_Count=9000/Yarn_Count
        return Yarn_Count
    
    elif sys1=='nm' and sys2=='ne':
        Yarn_Count=0.59*Yarn_Count
        return Yarn_Count
    
    elif sys1=='nm' and sys2=='nw':
        Yarn_Count=0.89*Yarn_Count
        return Yarn_Count
    
    elif sys1=='nw' and sys2=='tex':
        Yarn_Count=885.8/Yarn_Count
        return Yarn_Count
    
    elif sys1=='nw' and sys2=='dtex':
        Yarn_Count=8858/Yarn_Count
        return Yarn_Count
    
    elif sys1=='nw' and sys2=='den':
        Yarn_Count=7972/Yarn_Count
        return Yarn_Count
    
    elif sys1=='nw' and sys2=='nm':
        Yarn_Count=1.129*Yarn_Count
        return Yarn_Count
    
    elif sys1=='nw' and sys2=='ne':
       Yarn_Count=(2/3)*Yarn_Count
       return Yarn_Count 
    
    else:
        
        print("Your inputs are invalid!")



#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#===-================================
#------------data analysis--



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



