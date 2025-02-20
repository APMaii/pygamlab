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



