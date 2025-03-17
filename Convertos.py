
'''
Convertos.py :










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




def Convert_Coulomb_to_Electron_volt( coulomb):
    electron_volt = coulomb * 6.24e18 
    return electron_volt

def Convert_Electron_volt_to_Coulomb( electron_volt):
    coulomb = electron_volt / 6.24e18
    return coulomb




def Celsius_To_Kelvin (C):
    K = float(C + 273.15) 
    return K 


#vice versa_______________

def Kelvin_To_Celsius (K): 
    C = float(K - 273.15)
    return C
    


def convert_cm2_cs_to_m2_s (cm2_cs) :
    m2_s = float(cm2_cs / 100)  # Conversion factor for area: 1 m^2 = 10000 cm^2  # Conversion factor for time: 1 s = 100 cs  # (10**-2) Combined conversion factor for square meters per second to square centimeters per centisecond
    return m2_s


#vice versa_________________________

def convert_m2_s_to_cm2_cs (m2_s) :
    cm2_cs = float(m2_s * 100)
    return cm2_cs

    


def Atmosphere_to_Pascal (atm):
    Pa= float(atm * 101325)
    return Pa 

    
#vice versa___________________

def Pascal_to_Atmosphere (Pa):
    atm = float(Pa / 101325)
    return atm 




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






#convertor1
def Radians_to_degrees (num):
    '''
    
    This function is used for convert radians to degree
    '''
    degree=num*180/pi
    return degree

#convertor2
def Weightpercent_to_ppm (num):
    '''
    
    This function is used for convert weight percent to ppm
    '''
    ppm= num*10000
    return ppm


    

def M_to_mm (M):
    '''
    
    This function is usef for convert Meter to Milimeter
    '''
    mm=M*1000
    return mm



#15000000

def mm_to_M (mm):
    '''
    
    This function is usef for convert Meter to Milimeter
    '''
    M=mm/1000
    return M

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

#-----------------------------

def C_to_Fahrenheit(C):
    '''
    This function is usef for Centigrade_to_Fahrenheit
    '''
    F=C*1.8+32 #cantigerad ra zarbedar1.8 bealave32 mikonim ta  Fahrenheit be dast ayad
    return F #dar in ja az return estefademikonim ta khoroji bedast ayad

#---------------------------




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

#------------------------------------------------------------------------------
#convertor1
def sec_to_hour(t):
    t=t/3600
    return t

#convertor 1 reversed
def hour_to_sec(t):
    t=t*3600
    return t

#------------------------------------------------------------------------------
#convertor 2 
def Cm_to_um(x):
    a=x*10000
    return a

#convertor 2 reverse
def um_to_Cm(x):
    a=x/10000
    return a





#===============================
#1_Constants:







#===========================================
#2_Functions:







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






#convertors
def Meter_to_nanometer(value):
    return value * (10**9)
    
def Nanometer_to_meter(value):
    return value * (10 ** (-9))

def Centigrade_to_kelvin(value):
    return value + 237






#part3

def Flux_convertor1(a):
    '''
This function is used for convert ml/cm2.min to L/m2hr
    a : flow in ml/cm2.s unit
    Returns flow in L/m2hr unit
    '''
    b=a*600
    return b


def Flux_convertor2(a):
    '''
This function is used for convert L/m2hr to ml/cm2.min
    a : flow in L/m2hr unit
    Returns flow in ml/cm2.min unit
    '''
    b=a/600
    return b


def Pressure_convertor1(a):
    '''
This function is used for convert mmHg to bar
    a : Pressure in mmHg unit
    Returns Pressure in bar unit
    '''
    b=a/760
    return b


def Pressure_convertor2(a):
    '''
This function is used for convert bar to mmHg
    a : Pressure in bar unit
    Returns Pressure in mmHg unit
    '''
    b=a*760
    return b


def Temperature_convertor1(a):
    '''
This function is used for convert Celsius to Fahrenheit
    a : Temperature in Celsius unit
    Returns Temperature in Fahrenheit unit
    '''
    b=(a*1.8)+32
    return b

def Temperature_convertor2(a):
    '''
This function is used for convert Fahrenheit to Celsius
    a : Temperature in Fahrenheit unit
    Returns Temperature in Celsius unit
    '''
    b=(a-32)/(1.8)
    return b


def Velocity_convertor1(a):
    '''
This function is used for convert m/s to km/hr
    a : Velocity in m/s unit
    Returns Velocity in km/hr unit
    '''
    b=a*3.6
    return b


def Velocity_convertor2(a):
    '''
This function is used for convert km/hr to m/s
    a : Velocity in km/hr unit
    Returns Velocity in m/s unit
    '''
    b=a/3.6
    return b



#part2



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


###---------------Part2-------------------------------------------------------
##part2-section1*********
##polymer burning rate calculation function


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

    






#---------------------
#part2:

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




#====================================================================================================================================================

#Part 2

#1 related to additive manufacturing






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
    
    




#2(3 ta tabe...)
#__________________________________________________________________________________________




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



#------------------2nd--------------------------------------------------------

# 1st function calculates Latice parameter of crystal structure.

# FCC(face centered cubic) 
# BCC(body centered cubic) 
# SC(simple cubic) 
# HCP(hexagonal close pack)
# DC(diamond cubic)



#------------------------------------------------------------------------------   

# 2nd function calculates fracture toughness 
# based on applied stress and crack length and location


#------------------------------------------------------------------------------       

# 3rd function calculates wear rate.


#------------------------------------------------------------------------------

# 4th function calculates Vickers hardness.



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









#2:
    #2-1:






  
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


#------------------------------
#Part 2


    
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


import math


#Converters
def Fahrenheit_to_Kelvin(F):
    return (F - 32) * 5 / 9 + 273.15

def Kelvin_to_Fahrenheit(K):
    return (K - 273.15) * 9 / 5 + 32

def Meter_to_Angstrom(m):
    return m*1e10

def Angstrom_to_Meter(A):
    return A*1e-10

def Milimeter_to_Angstrom(mm):
    return mm*1e7

def Angstrom_to_Milimeter(A):
    return A*1e-7

def Nanometer_to_Angstrom(nm):
    return nm*10

def Angstrom_to_Nanometer(A):
    return A/10

def Micrometer_to_Angstrom(um):
    return um*10000

def Angstrom_to_Micrometer(A):
    return A/10000


#section three
def Kilometer_Per_Hour_To_Meter_Per_Second(kph):
    '''
    Parameters
    ----------
    kph: float
         number in kilometer per hour
    mps: float
         number in meter per second
    '''
    mps=kph*3.6
    return mps





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


def Kg_to_Ton(Kg):
        Ton=1000*Kg
        return Ton





def Cm_To_Inch(length = 1) :
    '''
    This function convert cm to inch.
    Convertor law :  inch = length/2.54
     
    Parameters
    ----------
    value : float, optional
        Convert input to inch. The default is 1.

    Returns
    -------
    return input that converted to inch.

    '''
    inch = length/2.54
    return inch
    

def Inch_To_Cm(length = 1) :
    '''
    This function convert inch to cm.
    Convertor law :  inch = length*2.54

    Parameters
    ----------
    length : float, optional
        Convert input to cm. The default is 1.

    Returns
    -------
    return input that converted to cm..

    '''
    cm = length*2.54
    return cm




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

---------------------------------------------------------
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


#**************PART3**********CONVERTOR************************************8
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


#************************************************************
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




def Convert_Annual_To_Monthly_Loss(annual_loss):
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



print('*******************************************')
test=Convert_Annual_To_Monthly_Loss(34500000)
print (test)




def Convert_Persion_Hours_To_Persion_Days(hours, daily_working_hours):
    '''
    

    Parameters
    ----------
    hours : float
        total working hours.
    daily_working_hours : float
        Legal daily working hours.



    '''
    if not str(hours).isdigit() or not str(daily_working_hours).isdigit() :
        print ('error! bad parameter!')
    return hours/daily_working_hours



def Molarity_to_Normality(Molarity,n):
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
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Normality_to_Molarity(Normality,n):
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
    

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Mass_to_Mole(Mass,Molar_Mass):
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def Mole_to_Mass(Mole,Molar_Mass):
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


#----------3----------
def Kg_to_Lbm(Kg):
    Lbm=0.4535*Kg
    return Lbm
   

def Lbm_to_Kg(Lbm):
    Kg=2.20462*Lbm
    return Kg


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
#            ---------------------------------------------------
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
 
#==============================================================================
def Decimal_To_Binary(Num_dec):
    Bin=0
    i=0
    while Num_dec!=0:
        r=Num_dec%2
        Bin=Bin+(r*(10**i))
        Num_dec=Num_dec//2
        i=i+1
    return Bin



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




def Centimeter_per_Minute_to_Meter_per_Hour_Welding_Speed_Converter(Centimeter_per_Minute):
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


def Meter_per_Hour_to_Centimeter_per_Minute_Welding_Speed_Converter(Meter_per_Hour):
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


# 3. Liter/Minute to CC/Second Welding Gas Flow Rate Converter
    
def Liter_per_Minute_to_CC_per_Second_Welding_Gas_Flow_Rate_Converter(Liter_per_Minute):
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


# 4.  CC/Second to Liter/Minute Welding Gas Flow Rate Converter
    
def CC_per_Second_to_Liter_per_Minute_Welding_Gas_Flow_Rate_Converter(CC_per_Second):
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


def Mm_year_to_Mils_year(milpy):
    """
    1mm/yr=39.37mpy
    Crossion rate
     """
    mpy=39.37*milpy
    return mpy
# Mm_year_to_Mils_year(int(input()))

def Mils_year_to_Mm_year(mpy):
    """
      1mm/yr=39.37mpy
      Crossion rate
    """
    Mm_year=mpy/39.37
    return Mm_year
# print(Mils_year_to_Mm_year(float(input())))



def Rockwell_to_Brinell(hrb):
    '''
    convert Rockwell hardness (HRB) to Brinell hardness (HB).
    
    Parameters
    
    hrb : float
        hardness in Rochwell scale.

    Returns float: Hardness in Brinell scale.
    

    '''
    hb = (hrb * 5.0) + 50
    return hb


    
def Brinell_to_Rockwell(hb):
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



def Horsepower_to_Watt (Horsepower):
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



def Watt_to_Horsepower (Watt) :
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




def Force_CGS_to_SI (Force_in_CGS):
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

def Force_SI_to_CGS (Force_in_SI) :
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


def Current_density_to_mpy(Current_density,density,masschange,valency):
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

def  Mpy_to_current_density(mpy,density,masschange,valency):
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




