
'''
Convertors.py :

This module provides converter functions for transforming values between different units of measurement.

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



import math

def Atmosphere_To_Pascal(atm):
    """
    Convert the given value in atmospheres to pascals.

    Parameters
    ----------
    atm : float
        The value in atmospheres to be converted.

    Returns
    -------
    float
        The equivalent value in pascals.
    """
    Pa = float(atm * 101325)
    return Pa


def Angstrom_To_Meter(A):
    """
    Convert the given value in angstroms to meters.

    Parameters
    ----------
    A : float
        The value in angstroms to be converted.

    Returns
    -------
    float
        The equivalent value in meters.
    """
    return A * 1e-10


def Angstrom_To_Millimeter(A):
    """
    Convert the given value in angstroms to millimeters.

    Parameters
    ----------
    A : float
        The value in angstroms to be converted.

    Returns
    -------
    float
        The equivalent value in millimeters.
    """
    return A * 1e-7


def Angstrom_To_Nanometer(Angstrom_value):
    
    '''
    This function converts Angstroms to Nanometers.

    Parameters
    ----------
    Angstrom_value: int or Float
        Value in angstroms (Å).
    
    Returns
    -------
    Nanometer_value: int or Float
        Equivalent value in Nanometers (nm).

    '''
    Nanometer_value= Angstrom_value/10
    return Nanometer_value 




def Angstrom_To_Micrometer(A):
    """
    Convert the given value in angstroms to micrometers.

    Parameters
    ----------
    A : float
        The value in angstroms to be converted.

    Returns
    -------
    float
        The equivalent value in micrometers.
    """
    return A / 10000


def Annual_To_Monthly_Loss(annual_loss):
    '''
    

    Parameters
    ----------
    annual_loss : int
        the annual loss of an Economic enterprise.

    Returns
    -------

        the monthly loss of an Economic enterprise.

    '''
    if not str(annual_loss).isdigit(): 
        print ('error! bad parameter!')
    return int(annual_loss/12)





def Byte_To_Kilobyte(b):
    """
    Convert the given value in bytes to kilobytes.

    Parameters
    ----------
    b : int
        The value in bytes to be converted.

    Returns
    -------
    float
        The equivalent value in kilobytes.
    """
    kb = 0.0009765625 * b
    return kb



def Bar_To_Pascal(bar):
    """
    Converts bar to pascal.

    Parameters
    ----------
    bar : float
        bar.

    Returns
    -------
    Pa : float
        pascal

    """
    
    Pa=bar*(10**(-5))
    return Pa

def Brinell_To_Rockwell(hb):
    '''
    convert Brinell hardness (HB) to Rockwell hardness (HRB)

    Parameters
    ----------
    hb : float
        hardness in Brinell scale.

    Returns float: Hardness in Rochwell scale.
   

    '''
    
    hrb = (hb - 50) / 5.0
    return hrb




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
    


def Celcius_To_Kelvin(Celcius):
    """
    Convert the given temperature in Celsius to Kelvin.

    Parameters
    ----------
    Celcius : float
        The temperature in Celsius to be converted.

    Returns
    -------
    float
        The equivalent temperature in Kelvin.

    Notes
    -----
    The temperature in Celsius is offset by 273.15 to convert to Kelvin.
    """
    Kelvin = Celcius + 273.15
    return Kelvin


def Coulomb_To_Electron_volt(coulomb):
    """
    Convert the given value in coulombs to electron volts.

    Parameters
    ----------
    coulomb : float
        The value in coulombs to be converted.

    Returns
    -------
    float
        The equivalent value in electron volts.

    Notes
    -----
    1 coulomb = 6.24e18 electron volts.
    """
    electron_volt = coulomb * 6.24e18
    return electron_volt


def Centigrade_To_Fahrenheit(C):
    """
    Convert the given temperature in Celsius (Centigrade) to Fahrenheit.

    Parameters
    ----------
    C : float
        The temperature in Celsius to be converted.

    Returns
    -------
    float
        The equivalent temperature in Fahrenheit.
    """
    F = C * 1.8 + 32
    return F


def CmHg_To_Pascal(P1):
    """
    Convert the given pressure in centimeters of mercury to pascals.

    Parameters
    ----------
    P1 : float
        The pressure in centimeters of mercury to be converted.

    Returns
    -------
    float
        The equivalent pressure in pascals.
    """
    P2 = P1 * 1333.22
    return P2


def Calories_To_Joules(cal):
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


def Centimeter_per_Minute_To_Meter_per_Hour_Welding_Speed_Converter(Centimeter_per_Minute):
    '''
    This function converts the Welding Speed from Centimeter per Minute to Meter per Hour.

    Parameters
    ----------
    Centimeter_per_Minute : float
        Centimeter_per_Minute is a unit for welding speed.

    Returns
    -------
    Meter_per_Hour is a unit for welding speed.

    '''     
 
    Meter_per_Hour=Centimeter_per_Minute/1.7
    return Meter_per_Hour


def Current_density_To_mpy(Current_density,density,masschange,valency):
    """
    

    Parameters
    ----------
    Current_density : float
        Current density .(microA/cm2)
    density : float
       Material Density (g/cm3).
    masschange : float 
        amount of matter already corroded (g)
    valency : intiger
       How positive is the charge of the Material

    Returns
    -------
   corrosion rate in mpy
   

    """
    corrosion_rate_mpy=Current_density*1e-6*31536000*(1/density)*masschange*400*(1/(valency*96500))
    return corrosion_rate_mpy

def CC_per_Second_To_Liter_per_Minute_Welding_Gas_Flow_Rate_Converter(CC_per_Second):
    '''
    This function converts the Welding Gas Flow Rate from CC per Second to Liter per Minute.

    Parameters
    ----------
    CC_per_Second : float
        CC_per_Second is a unit for gas flow rate in welding.

    Returns
    -------
    Liter_per_Minute is a unit for gas flow rate in welding.

    '''     
 
    Liter_per_Minute=CC_per_Second/16.67
    return Liter_per_Minute



def Degree_To_Radian(deg):
    """
    Converts values of angle from degree to radian.

    Parameters
    ----------
    deg : float
        The angle value in degrees.

    Returns
    -------
    float
        The angle value in radians.
    """
    rad = deg * 3.141592653589793 / 180
    return rad


def Decimal_To_Binary(Num_dec):
    """
    Converts a decimal number to its binary representation.

    Parameters
    ----------
    Num_dec : int
        The decimal number to convert.

    Returns
    -------
    int
        The binary representation of the decimal number.
    """
    Bin = 0
    i = 0
    while Num_dec != 0:
        r = Num_dec % 2
        Bin = Bin + (r * (10**i))
        Num_dec = Num_dec // 2
        i = i + 1
    return Bin


def Electronvolt_To_Joule(e_v):
    """
    Converts energy value from electronvolts to joules.

    Parameters
    ----------
    e_v : float
        The energy value in electronvolts.

    Returns
    -------
    float
        The energy value in joules.
    """
    Joule = e_v * 1.6022e-19
    return Joule


def Electron_volt_To_Coulomb(electron_volt):
    """
    Converts energy value from electronvolts to coulombs (assuming it refers to charge equivalent).
    Note: This conversion is based on the elementary charge.

    Parameters
    ----------
    electron_volt : float
        The energy value in electronvolts.

    Returns
    -------
    float
        The equivalent charge value in coulombs.
    """
    coulomb = electron_volt * 1.602e-19
    return coulomb


def Foot_To_Mile(ft):
    """
    Converts a length value from feet to miles.

    Parameters
    ----------
    ft : float
        The length value in feet.

    Returns
    -------
    float
        The length value in miles.
    """
    mi = 0.000189393939 * ft
    return mi


def Fahrenheit_To_Centigrade(F):
    """
    Converts a temperature value from Fahrenheit to Centigrade (Celsius).

    Parameters
    ----------
    F : float
        The temperature value in Fahrenheit.

    Returns
    -------
    float
        The temperature value in Centigrade (Celsius).
    """
    C = (F - 32) * 5/9
    return C





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

    Newton_Meters=Foot_Pounds*1.3558
    return Newton_Meters


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




def Force_CGS_To_SI (Force_in_CGS):
    '''
    

    Parameters
    ----------
    Force_In_CGS : float
        give your force value in CGS system.

    Returns
    -------
    SI : float
        return your force value in SI system.

    '''
    
    SI = "{:e}".format(Force_in_CGS * 1e-5)
    return SI

def Force_SI_To_CGS (Force_in_SI) :
    '''
    

    Parameters
    ----------
    Force_in_SI : float
        give your force value in SI system.

    Returns
    -------
    CGS : float
        return your force value in CGS system.

    '''
    
    CGS = "{:e}".format(Force_in_SI * 1e+5)
    return CGS





def Gram_To_Mole(g,MW):
    '''
    This function calaculates the eqivalent amount of substance of a compound  in mole(s) base on mass in gram(s).

    Parameters
    ----------
    g : float
        g is the mass of a compound in gram(s).
    MW : float
        MW is the Molecular weight of a compound (gram/mol).

    Returns
    -------
    Mole : float
        Mole is the eqivalent amount of substance of a compound in mole(s).

    '''
    Mole = g / MW
    return Mole



def Hour_To_Sec(t):
    """
    Converts a time value from hours to seconds.

    Parameters
    ----------
    t : float
        The time value in hours.

    Returns
    -------
    float
        The time value in seconds.
    """
    t = t * 3600
    return t

def Hertz_To_Rpm(a,/):
    '''
    A converter machine to convert frequency in Hertz(Hz) to frequency in rpm.
    Parameters
    ----------
    a : int or float
        frequency, Hertz(Hz).

    Returns
    b : int or float 
    frequency, revolution per minute (rpm)
    '''
    b=a*60
    return b


def Horsepower_To_Watt (Horsepower):
    '''
    

    Parameters
    ----------
    Horsepower : float
        give number in horsepower.

    Returns
    -------
    watt : float
        return your number in watt.

    '''
    Watt = "{:e}".format(Horsepower * 745.7)
    return Watt




def Inch_To_Centimeter(Inch):
    '''
    Parameters
    ----------
    Inch : float or int
        ne inch is equal to 2.54 centimeters.
        number per Inch unit.

    Returns
    -------
    Centimeter : float
        number per Centimeter unit.

    '''
    Centimeter=2.54*Inch
    return Centimeter

def Inch_To_Meter(In):
    """
    Converts a length value from inches to meters.

    Parameters
    ----------
    In : float
        The length value in inches.

    Returns
    -------
    float
        The length value in meters.
    """
    m = In / 39.3701
    return m



def Joules_To_Calories(J):
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

def Joules_Per_Minute_To_Kilowatt(Joules_Per_Minute):
    '''

    Parameters
    ----------
    Joules_Per_Minute : float
        number per Joules unit.

    Returns
    -------
    Kilowatt : float
        number per Kilowatt unit.

    '''
    Kilowatt=(Joules_Per_Minute)/60000
    return Kilowatt




def Kilogram_To_Pound(number_in_kilogram):
    '''
    This function converts the desired number from kilograms to pounds.

    Parameters
    ----------
    number_in_kilogram : int
        Number per kilogram.

    Returns
    -------
    pound : int
        Number per pound.

    '''
    pound=number_in_kilogram*2.2046
    return pound


def Kelvin_to_Celcius(Kelvin):
    """
    This function is used to convert Kelvin to Celsius.
    The temperature in Celsius is different from the temperature in Kelvin by 273.15.

    Parameters
    ----------
    Kelvin : float
        The temperature value in Kelvin.

    Returns
    -------
    float
        The temperature value in Celsius.
    """
    Celcius = Kelvin - 273.15
    return Celcius

def Kilogeram_Per_Cubic_Meter_To_Pounds_Per_Cubic_Inch(KgPerCubicMeter):
    """
    Converts a density value from kilograms per cubic meter to pounds per cubic inch.

    Parameters
    ----------
    KgPerCubicMeter : float
        The density value in kilograms per cubic meter.

    Returns
    -------
    float
        The density value in pounds per cubic inch.
    """
    L = KgPerCubicMeter * 0.0000361273
    return L


def KiloMeter_To_LightYear(km):
    """
    Converts a distance value from kilometers to light-years.

    Parameters
    ----------
    km : float
        The distance value in kilometers.

    Returns
    -------
    float
        The distance value in light-years.
    """
    ly = km / 9460730472801.1
    return ly



def Kmph_To_Mps(V1):
    """
    This function is used to convert kilometers per hour to meters per second.

    Parameters
    ----------
    V1 : float
        The speed value in kilometers per hour.

    Returns
    -------
    float
        The speed value in meters per second.
    """
    V2 = V1 / 3.6
    return V2


def Kilobyte_To_Byte(kb):
    """
    Converts a data size value from kilobytes to bytes.

    Parameters
    ----------
    kb : float
        The data size value in kilobytes.

    Returns
    -------
    float
        The data size value in bytes.
    """
    b = 1024 * kb
    return b



def Kilometer_Per_Hour_To_Meter_Per_Second(kph):
    """
    Converts a speed value from kilometers per hour to meters per second.

    Parameters
    ----------
    kph : float
        The speed value in kilometers per hour.

    Returns
    -------
    float
        The speed value in meters per second.
    """
    mps = kph / 3.6
    return mps


def Kg_To_Ton(Kg):
    """
    Converts a mass value from kilograms to metric tons.

    Parameters
    ----------
    Kg : float
        The mass value in kilograms.

    Returns
    -------
    float
        The mass value in metric tons.
    """
    Ton = Kg / 1000
    return Ton


def Kg_To_Lbm(Kg):
    """
    Converts a mass value from kilograms to pounds (lbm).

    Parameters
    ----------
    Kg : float
        The mass value in kilograms.

    Returns
    -------
    float
        The mass value in pounds (lbm).
    """
    Lbm = Kg * 2.20462
    return Lbm




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



def LightYear_To_KiloMeter(ly):
    """
    Converts a distance value from light-years to kilometers.

    Parameters
    ----------
    ly : float
        The distance value in light-years.

    Returns
    -------
    float
        The distance value in kilometers.
    """
    km = ly * 9460730472801.1
    return km


def Lbm_To_Kg(Lbm):
    """
    Converts a mass value from pounds (lbm) to kilograms.

    Parameters
    ----------
    Lbm : float
        The mass value in pounds (lbm).

    Returns
    -------
    float
        The mass value in kilograms.
    """
    Kg = Lbm / 2.20462
    return Kg



def Liter_per_Minute_To_CC_per_Second_Welding_Gas_Flow_Rate_Converter(Liter_per_Minute):
    '''
    This function converts the Welding Gas Flow Rate from Liter per Minute to CC per Second.

    Parameters
    ----------
    Liter_per_Minute : float
        Liter_per_Minute is a unit for gas flow rate in welding.

    Returns
    -------
    CC_per_Second is a unit for gas flow rate in welding.

    '''     
 
    CC_per_Second=Liter_per_Minute*16.67
    return CC_per_Second




#--------m------


def Meter_To_MilliMeter(meter):
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



def MilliMeter_To_Meter (milimeter):
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



def Micrometer_To_Nanometer(micrometer):
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
    Nanometer=micrometer*1000
    return Nanometer

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
    
    Pascal=Megapascal/1000000
    return Pascal


def Mps_To_Kmph(V1):
    """
    This function is used to convert meter per second to kilometer per hour.

    Parameters
    ----------
    V1 : float
        The speed value in meters per second.

    Returns
    -------
    float
        The speed value in kilometers per hour.
    """
    V2 = V1 * 3.6
    return V2


def Mile_To_Foot(mi):
    """
    Converts a length value from miles to feet.

    Parameters
    ----------
    mi : float
        The length value in miles.

    Returns
    -------
    float
        The length value in feet.
    """
    ft = 5280 * mi
    return ft


def Meter_To_inch(m):
    """
    Converts a length value from meters to inches.

    Parameters
    ----------
    m : float
        The length value in meters.

    Returns
    -------
    float
        The length value in inches.
    """
    In = m * 39.3701
    return In



def Miller_To_Millerbrove(u,v,w):
    
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
  

def Millerbrove_To_Miller(l,m,n,o):

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


def Meter_To_Angstrom(m):
    """
    Converts a length value from meters to Angstroms.

    Parameters
    ----------
    m : float
        The length value in meters.

    Returns
    -------
    float
        The length value in Angstroms.
    """
    return m * 1e10


def Milimeter_To_Angstrom(mm):
    """
    Converts a length value from millimeters to Angstroms.

    Parameters
    ----------
    mm : float
        The length value in millimeters.

    Returns
    -------
    float
        The length value in Angstroms.
    """
    return mm * 1e7


def Micrometer_To_Angstrom(um):
    """
    Converts a length value from micrometers to Angstroms.

    Parameters
    ----------
    um : float
        The length value in micrometers.

    Returns
    -------
    float
        The length value in Angstroms.
    """
    return um * 10000




def Meter_Per_Second_To_Kilometer_Per_Hour(mps):
    '''
    Parameters
    ----------
    mps: float
         number in meter per second
    kph: float
         number in kilometer per hour
    '''
    kph=mps/3.6
    return kph






def Mole_To_Gram(mol,MW):
    '''
    This function calaculates the eqivalent mass of a compound in gram(s) base on amount of substance in mole(s).

    Parameters
    ----------
    mol : float
        mol is the eqivalent amount of substance of a compound in mole(s).
    MW : float
        MW is the Molecular weight of a compound (gram/mole).

    Returns
    -------
    g : float
        g is the eqivalent mass in of a compound in in gram(s).

    '''
    g = mol * MW
    return g





def Mass_To_Mole(Mass,Molar_Mass):
    '''
    

    Parameters
    ----------
    Mass : float
        The mass of substance(g).
    Molar_Mass : float
        The mass of one mole of substance (g/mol).

    Returns
    -------
    Mole: int

    '''
    Mole=Mass/Molar_Mass
    return(Mole)



def Mole_To_Mass(Mole,Molar_Mass):
    '''
    

    Parameters
    ----------
    Mole : int
        
    Molar_Mass : float
        The mass of one mole of substance (g/mol).

    Returns
    -------
    Mass (g) : Float.

    '''
    Mass=Mole*Molar_Mass
    return(Mass)


def Mpa_To_Psi(Num_Mpa,/):
    '''
    

    Parameters
    ----------
    
    Num_Mpa : float
        Megapascals=Newton per square millimetre

    Returns
    -------
    Psi : float
        Psi=Pounds force per square inch 

    '''
    Psi=Num_Mpa*145
    return Psi


def Meter_per_Hour_To_Centimeter_per_Minute_Welding_Speed_Converter(Meter_per_Hour):
    '''
    This function converts the Welding Speed from Meter per Hour to Centimeter per Minute.

    Parameters
    ----------
    Meter_per_Hour : float
        Meter_per_Hour is a unit for welding speed.

    Returns
    -------
    Centimeter_per_Minute is a unit for welding speed.

    '''     
 
    Centimeter_per_Minute=Meter_per_Hour*1.7
    return Centimeter_per_Minute



def Mm_year_To_Mils_year(milpy):
    """
    Converts a corrosion rate from millimeters per year (mm/yr) to mils per year (mpy).
    1 mm/yr = 39.37 mpy

    Parameters
    ----------
    milpy : float
        The corrosion rate in millimeters per year.

    Returns
    -------
    float
        The corrosion rate in mils per year.
    """
    mpy = 39.37 * milpy
    return mpy


def Mils_year_To_Mm_year(mpy):
    """
    Converts a corrosion rate from mils per year (mpy) to millimeters per year (mm/yr).
    1 mm/yr = 39.37 mpy

    Parameters
    ----------
    mpy : float
        The corrosion rate in mils per year.

    Returns
    -------
    float
        The corrosion rate in millimeters per year.
    """
    Mm_year = mpy / 39.37
    return Mm_year



def  Mpy_To_current_density(mpy,density,masschange,valency):
    """
    

    Parameters
    ----------
    mpy : float
        corrosion rate in mpy
    density : float
        materails density 
    masschange : float
        amount of mass corroded 
    valency : int
        how positive is the charge

    Returns
    -------
    Current density 

    """
    Current_density=(mpy*1e6*density*2.5*valency*96500)/(31536000*masschange*1000)
    return Current_density



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





def Nanometer_To_Micrometer(nanometer):
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
    Micrometer=nanometer/1000
    return Micrometer


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
 


def Normality_To_Molarity(Normality,n):
    '''
    

    Parameters
    ----------
    Normality : float
    n : int
        Number of moles.

    Returns
    -------
    Molarity.

    '''
    Molarity=Normality/n
    return(Molarity)
    


def Nanometer_To_Angstrom(Nanometer_value):
    
    '''
    This function converts Nanometers to Angstroms.
    1 Nanometer(nm)= 10 Angstroms(Å)

    Parameters
    ----------
    Nanometer_value: int or float
        Value in Nanometers(nm).
    
    Returns
    -------
    Angstrom_value: int or float
        Equivalent value in Angstroms(Å).

    '''
    Angstrom_value= Nanometer_value*10
    return Angstrom_value

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
    
    Foot_Pound=Newton_Meters*0.7376
    return Foot_Pound





def Pounds_Per_Cubic_Inch_To_Kilogeram_Per_Cubic_Meter(LbPerCubicInch):
    """
    Converts a density value from pounds per cubic inch to kilograms per cubic meter.

    Parameters
    ----------
    LbPerCubicInch : float
        The density value in pounds per cubic inch.

    Returns
    -------
    float
        The density value in kilograms per cubic meter.
    """
    Kg = LbPerCubicInch * 27679.9
    return Kg



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
    
    Megapascal=1000000*Pascal
    return Megapascal



  
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



def Pascal_To_Atmosphere(Pa):
    """
    Converts a pressure value from Pascals to atmospheres.

    Parameters
    ----------
    Pa : float
        The pressure value in Pascals.

    Returns
    -------
    float
        The pressure value in atmospheres.
    """
    atm = float(Pa / 101325)
    return atm


def Percentages_To_Moles(total, percentages):
    """
    Calculates the number of moles of each component in a mixture given their percentages by weight and the total weight.

    Parameters
    ----------
    total : float
        The total weight of the mixture.
    percentages : dict
        A dictionary where keys are the names of the materials and values are their weight percentages.

    Returns
    -------
    dict
        A dictionary where keys are the names of the materials and values are their corresponding number of moles.
    """
    # Define the molecular weight of the composite mixture
    molar_weight = {'TEGDMA': 156.27, 'BIS_GMA': 512.67, 'UDMA': 398.48,
                    'Silica dioxide': 60.08, 'Barium silicate': 233.39, 'Zirconium dioxide': 123.22
                    }

    # Calculate the moles of each material
    moles = {}
    for material, percent in percentages.items():
        moles[material] = (percent / 100) * (total / molar_weight[material])

    return moles



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








def Pascal_To_CmHg(P1):
    """
    This function is used to convert Pascal to centimeter mercury.

    Parameters
    ----------
    P1 : float
        The pressure value in Pascals.

    Returns
    -------
    float
        The pressure value in centimeters of mercury (cmHg).
    """
    P2 = P1 / 1333.22
    return P2


def Ppm_To_Percent(a):
    """
    Converts a concentration value from parts per million (ppm) to percent.

    Parameters
    ----------
    a : float
        The ion concentration in ppm in brine.

    Returns
    -------
    float
        The ion percent in brine.
    """
    b = a / 10000
    return b


def Percent_To_Ppm(a):
    """
    Converts a concentration value from percent to parts per million (ppm).

    Parameters
    ----------
    a : float
        The ion percent in brine.

    Returns
    -------
    float
        The ion concentration in ppm in brine.
    """
    b = a * 10000
    return b


def Pascal_To_Torr(pa):
    """
    This function converts Pascal to Torr.

    Parameters
    ----------
    pa : float
        The pressure value in Pascals.

    Returns
    -------
    float
        The pressure value in Torr.
    """
    torr = pa / 133.322
    return torr




def Pascal_To_Bar(Pa):
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


def Psi_To_Mpa(Num_Psi,/):
    '''
    

    Parameters
    ----------
    
    Num_Psi : float
        Psi = Pounds force per square inch 

    Returns
    -------
    Mpa : float
        Megapascals=Newton per square millimetre

    '''
    Mpa=Num_Psi*(1/145)
    return Mpa



 





def Pound_To_Kilogram(number_in_pound):
    '''
    This function converts the desired number from pounds to kilograms.

    Parameters
    ----------
    number_in_pound : int
        Number per pound.

    Returns
    -------
    kilogram : int
        Number per kilogram.

    '''
    kilogram=number_in_pound/2.2046
    return kilogram







def Molarity_To_Normality(Molarity,n):
    '''
    

    Parameters
    ----------
    Molarity : float
    n : int
        Number of moles.

    Returns
    -------
    Normality.

    '''
    Normality=Molarity*n
    return(Normality)





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


def Sec_To_Hour(t):
    """
    Converts a time value from seconds to hours.

    Parameters
    ----------
    t : float
        The time value in seconds.

    Returns
    -------
    float
        The time value in hours.
    """
    t = t / 3600
    return t


def Radians_To_Degrees(num):
    """
    This function is used to convert radians to degrees.

    Parameters
    ----------
    num : float
        The angle value in radians.

    Returns
    -------
    float
        The angle value in degrees.
    """
    degree = num * 180 / math.pi
    return degree

def Rpm_To_Hertz(b,/):
    '''
   A converter machine to convert frequency in rpm to frequency in Herta(Hz).
    Parameters
    ----------
    b : int or float
        frequency, revolution per minute (rpm).

    Returns
    a, frequency, Hertz(Hz)

    '''
    a=b/60
    return a



def Torr_To_Pascal(torr):
    """
    This function converts Torr to Pascal.

    Parameters
    ----------
    torr : float
        The pressure value in Torr.

    Returns
    -------
    float
        The pressure value in Pascals.
    """
    pa = torr * 133.322
    return pa


def Viscosity_To_Poise(pa_s):
    """
    Converts dynamic viscosity from Pascal-seconds (Pa·s) to Poise (P).

    Parameters
    ----------
    pa_s : float
        The dynamic viscosity in Pascal-seconds.

    Returns
    -------
    float
        The dynamic viscosity in Poise.
    """
    poise = pa_s * 10
    return poise


def Viscosity_To_Pas(poise):
    """
    Converts dynamic viscosity from Poise (P) to Pascal-seconds (Pa·s).

    Parameters
    ----------
    poise : float
        The dynamic viscosity in Poise.

    Returns
    -------
    float
        The dynamic viscosity in Pascal-seconds.
    """
    pa_s = poise / 10
    return pa_s

    






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




def Weightpercent_To_ppm(num):
    """
    This function is used to convert weight percent to ppm (parts per million).

    Parameters
    ----------
    num : float
        The concentration in weight percent.

    Returns
    -------
    float
        The concentration in ppm.
    """
    ppm = num * 10000
    return ppm


def Watt_To_Horsepower (Watt) :
    '''
    

    Parameters
    ----------
    Watt : float
        give number in Watt.

    Returns
    -------
    Horsepower : float
        return number in Horsepower.

    '''
    Horsepower = "{:e}".format(Watt / 745.7)
    return Horsepower






    




