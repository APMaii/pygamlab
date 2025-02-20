'''

contributers

A2-Mohammad-Shafiee
A2_Sam_Beigy
A2_Najmeh_ Azizizadeh
A2_Naghmeh_Gholamalizadeh2
A2_Mohamadali_Rabiei
A2_MINA_MOHAMMADI
A2_MARZIEH_ALIDADI
A2_Mahdi_Nematollahi
A2_Fatemeh-Babaei
A2_Elnaz_Askari
A2_Behnaz_Hadi
A2_ARSHAWN_REZVANIPOUR
A2_Amin_Abedi
A1-Mohammad-Shafiee
A1_Zahra_Tamimi
A1_pariya_isvandi.py
A1_Najmeh_ Azizizadeh.py
A1_Naghmeh_Gholamalizadeh.py
A1_Marziyeh_Hatami-2.py
A1_MARZIEH_ALIDADI-2
A1_Mahdi_Nematollahi
A1_Fatemeh Babaei
A1_Behnaz_Hadi-2
A1_ARSHAWN_REZVANIPOUR
A1_Amin_Abedi
A1_Ali_Moftakharzadeh


'''

import math
### constnats
me = 9.11e-31           # Electron Rest Mass Kg
mp = 1.67264e-27        # Proton Rest Mass kg 
e  = 1.60218e-19        # Elementary Charge C(coulomb)
c  = 2.99792e8          # Speed of Light in Vacuum m/s
h  = 6.62617e-34        # Planck Constant J-s(joule-seconds)
ħ  = 1.05458e-34        # Reduced Planck Constant (ħ = h/2π)
k  = 1.38066e-23        # Boltzmann Constant J/K (joules per kelvin)
eV = 1.60218e-19        # Electron Volt J (joules)

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

#Functions
def Thermal_Voltage(Temperature, Unit = 'K'):
    '''    
    Calulation of thermal voltage value
    Parameters
    ----------
    Temperature : 
        Temperature.
    Unit : The unit of temperature
        K or k: Kelvin
        F or f: fahrenheit
        C or c: celsius
        The default is 'K'.
    Returns
    -------
    Calulated thermal voltage(KT/q)
    '''
    if Unit.upper() == 'F':
       #Temperature = (Temperature - 32) * 5 / 9 + 273.15 
       Temperature = Fahrenheit_to_Kelvin(Temperature)
    elif Unit.upper() == 'C':
      Temperature = Temperature + 273.15
    #return 1.38066e-23*Temperature/1.60218e-19  
    return k*Temperature/e  


def Diod_Forward_Voltage(Temp, I, Is):
    '''
    Calculation of forward voltage of a diode
    Parameters
    ----------
    Temp : float
        diode temperature in celsius.
    I : float
        diode current(A).
    Is : float
        the reverse saturation current(A).

    Returns
    -------
    the voltage across the diode(V).
    '''
    Vt = Thermal_Voltage(Temp, 'C')
    Vd = Vt*math.log(I/Is+1)
    return Vd


#section one
E0=8.8541878128*10**(-12) #vacuum permittivity
u0=1.256637061436*10**(-6) #vacuum permeability
h=6.626068*10**(-34) #planck constant
k=1.380649*10**(-23) #boltzman constant
Ke=8.9879*10**(9) #coulomb constant
#======================================================================

#section two
def Power_Spectrum_Density(frequency,T):
    '''
    Parameters
    ----------
    frequency: float
               the frequency that power spectrum density is required at
    T: float
       the temperature of material in kelvin

    Returns
    -------
    s: float
       power spectrum density of thermal noise           
    '''
    s=(frequency*h)/2*(2.7**((h*frequency)/(k*T))-1)
    return s
    

def Voltage_Standing_Wave_Ratio(Vmax,Vmin):
    '''
    Parameters
    ----------
    Vmax: float
          the highest voltage measured in transition line
    Vmin: float
              the lowest voltage measured in transition line
    Returns
    -------
    Row: float
         its the domain of reflection coefficient of transition line
    '''
    s=Vmax/Vmin
    Row=(s-1)/(s+1)
    return Row


#=========================================================================


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



#part1=============================================================================================================================================
#Gravitational_force_constants_numbers

g=6.673*10**-11                         #(N.m**2/kg**2)Gravitational_constant
G=6.673*10**-11                         #(N.m**2/kg**2)Gravitational_constant
m_earth=5.972*(10**24)                  #(kg)Mass_of_Earth
m_sun=1.989*(10**30)                    #(kg)Mass_of_the_Sun
m_mars=6.41693*(10**23)                 #kg)Mass_of_the_Mars
r_earth_sun=1.496*(10**11)              #(m)Average_distance_from_Earth_to_the_Sun
r_mars_sun=220.14*(10**9)               #(m)Average_distance_from_Mars_to_the_Sun
F=G*((m_earth*m_sun)/r_earth_sun**2)    #(N)Force_of_gravity
F2=g*((m_mars*m_sun)/r_mars_sun**2)     #(N)Force_of_gravity2

#part2=============================================================================================================================================
#Gravitational_force_formula

def Gravitational_force(G,FT):
    '''
    Parameters
    ----------
    G : float
        The gravitational force is a force that attracts any two objects with mass.
    FT : float
        The force of gravity varies with latitude and increases from about 9.780 m/s2 at the Equator to about 9.832 m/s2 at the poles.

    Returns
    -------
    F : float
        The force of gravity is also called as Newton's law of gravitation. Mathematically, F = GMmr2.
        where F = force of gravity, G = gravitational constant, M = mass of one object, m = mass of other object, r = distance between two objects.

    '''
    if FT=='Force of gravity total':
        F=G*((5.972*(10**24)*1.989*(10**30))/((1.496*(10**11))**2))
        return F
    
def Gravitational_force_formula(g,m_mars,m_sun,r_mars_sun):
    '''
    

    Parameters
    ----------
    g : float
        The gravitational force is a force that attracts any two objects with mass.
    m_mars : float
        Mars' mass is 6.42 x 1023 kilograms, about 10 times less than Earth.
        This affects the force of gravity. Gravity on Mars is 38 percent of Earth's gravity,
        so a 100-pound person on Earth would weigh 38 pounds on Mars
    m_sun : float
        The sun has a mass of 1.9891x1030 kg = 4.384x1030 lb = 2.192x1027 tons, or a mass 333,000 times that of the Earth.
        The radius of the Sun is 696,265,000 meters = 696,265 km = 432,639 mi or a radius 109 times that of the Earth.
    r_mars_sun : float
        Mars is about 128 million miles (206 million km) from the sun,
        and at its farthest distance (aphelion) Mars is about 154 million miles (249 million km) from the sun.

    Returns
    -------
    F2 : float
        The force of gravity is also called as Newton's law of gravitation.
        Mathematically, F = GMmr2.
        where F = force of gravity, G = gravitational constant, M = mass of one object, m = mass of other object, r = distance between two objects.

    '''
    if g=='gravity':
       F2=g*((m_mars*m_sun)/r_mars_sun**2)
       return F2
    
#part3=============================================================================================================================================

    def Kg_to_Ton(Kg):
        Ton=1000*Kg
        return Ton
    
    
    
    
    Kg_to_Ton(1)


    


#%% Part1: Introducing five constants 

B  = 1.38e-23          # Boltzman constant  J/k
EC = 1.602e-19         # Electron charge  C
Er = 2.81792e-15       # Electron radius  m
FC = 9.648e4           # Faraday constant C/mol
PC = 6.626e-34         # Plank constant Js 
AG = 9.8               # Acceleration gravity  m/s^2


#%% Part2: Two functions 
#=============  Solving differential equation based Euler method ==============   
def Euler_Diff_Solver(a , b, xf , h , x0 =0 ,y0 =0 ):
    '''
    This function solve linear differential equation of the following form: 
        yprim = dy/dx = f(x,y) = a*x^2 + b*x
    by the Euler's method that is a numerical method for 
    approximating differential equations.    

    Parameters
    ----------
    a : float
        Polynomial coefficient..
    b : float
        Polynomial coefficient.
    xf : float
        Independent variable's final value.
    h : float
        Step size.        
    x0 : float, optional
        Independent variable's initial value. The default is 0.
    y0 : float, optional
        Initial function value at x0 (y0 = y(x0)). The default is 0.

    Raises
    ------
    TypeError
        The quadratic coefficient of the equation should not be zero.

    Returns
    -------
    y : float
        Returns function value at desired point (xf) or y(xf).

    '''
    
    if a == 0:
        raise TypeError('The coefficient a should not be zero')
        return
    varL =  xf - x0
    varLL = int(varL/h)
    for i in range(varLL):
        y = Euler_Method(a,b,x0,y0,h)
        y0 = y
        x0 = x0 +h
    return y   
   
    
def Euler_Method(a,b, h, x0,y0):
    
    '''
    Euler's method law: y = y0 + h *yprim 

    Parameters
    ----------
    a : float
        Polynomial coefficient.
    b : float
        polynomial coefficient.
    x0 : float
        Initial value.
    y0 = y(x0) : float
        function value.
    h : float
        Step size.

    Returns
    -------
    This function returns function value (y) at specific value (x).

    '''
    var1 = x0;
    yprim = a*(var1**2)+b*var1;
    y = y0 + h *yprim # Euler's method
    return y
   

#%% Part3: Length Convertor
#=====================  Convert cm to inch  ===================================

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



===========================
####ثابت بولتزمن 
######ثابت پلانک 
##ثابت جهانی گاز ها 
#####ثابت گرانش 
####عدد آوگادرو
k=1.380649*(10**-23)(J.K^-1)
h=6.62607015*(10**-34)(kg.m^2.S-1)
R=8.134(J.K^-1.mol^-1)
G=6.67*(10**-11)(m^3/Kg*S^2)
NA=6.02214076*(10**23)

===========================================
def Calculate_Pipe_Heat_Transfer(mu,muW,rho,Cp,Pr,K,u,d,l,Tw,Tb1):
    '''
    Parameters
    ----------
    mu : float
        Fluid dynamic viscosity (pa.s).
    muW : float
        Dynamic viscosity at wall temperature (pa.s).
    rho : float
        fluid density (kg/m³).
    Cp : float
        DESCRIPTION.
    Pr : float
        Prantel number.
    K : flout
        Thermal conductivity(W/m.K).
    u : float
       fluid velocity (m/s).
    d : float
        Pipe diameter (m).
    L : float
        Pipe length (m).
    Tw : float
        Pipe wall temperature.
    Tb1 : float
        Average fluid temperature().

    Returns
    Reynolds number, Nusselt number, convective heat transfer coefficient, total heat transfer, and outlet temperature.
    -------
    None.

    '''
########################
#توضیحات کد: ابتدا عدد رینولدز محاسبه می‌گرد، اگر که کمتر از 2000 باشد در این صورت جریان آرام است .
# برای جریان آرام داخل لوله اگر لوله طویل باشد عدد ناسلت 3.66 می شود و اکر لوله طویل نباشد از فرمول ناسلت نوشته شده (1) استفاده می گردد . 
#و اگر عدد رینولدز در بازه ی 2500-1.25ضرب در 10 به توان 5 باشد و عدد پرانتل در بازه ی ذکر شده باشد ،از فرمول 2 استفاده می شود در صورتی که برای n تگر دمای جداره بیشتر از دمای سیال باشد برابر 0.4 .
# حالا می خواهیم با استفاده از عدد ناسلت ضریب انتقال حرارت جابه جایی را بدست آوریم .
# و ضریب انتقال حرارنت جابه جایی را حساب نموده و با استفاده از تساوی آخر دمای خروجی سیال بدست می آید .
########################
    
    Re=(rho*u*d)/mu
    if Re<2000:
        print('flow type laminar')
        if (Re*Pr*(d/l))>10:
            Nu=1.86*(Re*Pr)**(1/3)*(d/l)**(1/3)*(mu/muW)**0.14#####برای لوله های طویل صادق نیست(1) .
        else:
            Nu=3.66####برای لوله های طویل و دمای ثابت در جداره لوله 
    elif 2500<Re<(1.25*(10**5)) and 1.5<Pr<100:
        print('flow type Turbulent')
        if Tw>Tb1:
            n=0.4
        else:
            n=0.3
        Nu=(0.023*(Re**0.8)*(Pr**n))##(2)
    else:
        return None
    h = (Nu * K) / d
    A=3.14*d*l
    Tb2 =Tb1+(h*A*(Tw-Tb1))/(rho*Cp*u) 
    Q=((rho*A*u)*Cp*(Tb2-Tb1))
    return {Re,Nu,h,Q,Tb2}   
        
results=Calculate_Pipe_Heat_Transfer(0.000417,0.00035,985,4.18,3.02,0.651,0.02,0.0254,3,80,60)



##############################################################
import math
def Heat_Exchanger_Transfer(U,Th1,Th2,Tc1,Tc2,C,dot_m):
    '''
    

    Parameters
    ----------
    U : float
        Overall heat transfer coefficient.
    Th1 : float
        Hot fluid inlet temperature.
    Th2 : float
        Hot fluid outlet temperature.
    Tc1 : float
        Cold fluid inlet temperature.
    Tc2 : float
        Cold fluid outlet temperature.
    C : float
        Special heat.
    dot_m : float
        Debbie Jeremy.
        

    Returns(delta_T_LMTD,Q,A)
    -------
    None.
    

    '''
    
    delta_T_LMTD=((Th1-Tc2)-(Th2-Tc1))/math.log((Th1-Tc2)/(Th2-Tc1))##Logarithmic mean temperature
    global Q
    Q=dot_m*C*(Tc2-Tc1)
    global A
    A=Q/(U*delta_T_LMTD)###The heat exchange surface used
    return delta_T_LMTD,A,Q


result= Heat_Exchanger_Transfer(320,110,75,35,75,4180,68)
print(Q,A)

===============================================================

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
energy=Joules_Per_Minute_To_Kilowatt(11371000)
print (energy)
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
length= Inch_To_Centimeter(20)
print(length)
#============================================== 

#===============================================

NA = 6.02214076 * (10 ** 23)     #Avogadro constant (1/mol)
h = 6.62607015 * (10 ** -34)     #Planck constant (J/Hz)
k = 1.380649 * (10 ** -23)       #Boltzmann constant (J/K)
F = 96485.3321233100184          #Faraday constant (C/mol)
c = 2.99792458 * (10 ** 8)       #Light speed in vacuum (m/s)

#===============================================

def Austenite_Martensite_VC(C):
    '''
    This function calaculates the volume change of a unit cell in Austenite to Marteniste transformation.

    Parameters
    ----------
    C : float
        C is the percentage of carbon in chemical composition of steel.

    Returns
    -------
    VC : float
        VC is the percentage of volume change of a unit cell in Austenite to Marteniste transformation.

    '''
    a0 = 3.548 + (0.044 * C)    #Austenite lattice parameter
    V0 = (a0 ** 3) / 2          #Volume of Austenite unit cell (FCC)
    a = 2.861 - (0.013 * C)     #Martensite lattice parameter
    c = 2.861 + (0.116 * C)     #Martensite lattice parameter
    V=c * (a ** 2)              #Volume of Martensite unit cell (BCT)
    Delta_V = V - V0 
    VC = (Delta_V / V0) * 100
    return VC



def Standard_Deviation(a):
    '''
    This function calculates the standard deviation of a list of numbers.

    Parameters
    ----------
    a : list
        a is a list of numbers.

    Returns
    -------
    SD : float
        SD is the standard deviation of a list of numbers.

    '''
    E = 0
    F = 0
    
    for i in a:
        E = E + i       #sum of numbers
    M = E / len(a)      #Mean
    
    for i in a:
        D = (i - M) ** 2
        F = F + D
    SD2 = (F / (len(a) - 1))
    SD = math.sqrt(SD2)
    return SD

#===============================================

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


import math
#***********PART1**********CONSTANT NUMBERS*******************************
N=6.02214*math.pow(10,23)
#Avogadro Number, 1/mole
R=8.3145
#Universal Gas Constant, j/(mol.K)
Kb=1.380649*math.pow(10,-23)
#Boltzmann Constant, j/K
Ksb=5.67*math.pow(10,-8)
#Stephan-Boltzmann Constant, W/(k^4.m^2)
h=6.6236*math.pow(10,-36)
#Planck Constant, j.s


#***********PART2***********FUNCTIONS********************
#William, Landel, Ferrry (WLF)
def WLF(T,Tg,/):
    '''
    The WLF equation is a procedure for shifting data for amorphous polymers obtained at elevated temperatures to a reference temperature. 

    Parameters
    ----------
    T : int or float
        Temperature, K or degree celsius, Tg<T<Tg+100..
    Tg : int or float
        Glass Transition Temperature, K or degree celsius.
   
    Returns
    -------
    aT : int or float
    shift factor.

    '''
    b=T-Tg 
    c=-17.44*b
    d=51.6+b
    e=c/d
    aT=math.pow(10,e)
    return aT
        

#Cohen Equation
def Cohen(m,r,/):
    '''
  Cohen equation is used to predict the yield stress of a polymeric blend  containing a rubbery dispersion phase. 

    Parameters
    ----------
    m : int or float
        The yield stress of polymeric matrix, N/m^2 or Pa or ...
    r : float
        Volume fraction of rubbery phase, 0<r<1.

    Returns
    -------
    b : int or float
    The yield stress of polymeric blend, N/m^2 or Pa or ...

    '''
    a=(1-1.21*math.pow(r,(2/3)))
    b=m*a
    return b
    
        
#Critical diameter of rubbery particles
def Critical_Diameter(d,r,/):
    '''
    This equation predicts the critical diameter of rubber particles toughening a polymeric matrix.
    Parameters
    ----------
    d : int or float
        critical distance between rubbery particles, angstrom or mm or nm or .....
    r : float
        Volume fraction of rubbery phase, 0<r<1.

    Returns
    -------
    dc : int or float
    the critical diameter of rubber particles

    '''
    a=6*math.pow(r,(1/3))
    b=a-1
    c=3.14/b
    dc=d/c
    return dc


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





#---------------------
#constants
#---------------------
LEARNING_RATE = 0.01
MAX_DEPTH = 10
DEFAULT_TIMEOUT = 30
CACHE_SIZE = 1024



#part 2
'''
implementing two functions in my career, software engineering and enterprise architecture

'''

#------------------------------
#first function: 
#insertion sort 
#------------------------------

def Insertion_Sort (non_sorted_list):
    
    # key = non_sorted_list[0]
    for i in range(1, len(non_sorted_list)):
        key = non_sorted_list[i]
        # print ('key:',key)
        # print ('i:',i)
        j=i-1
        while j>=0 and key<non_sorted_list[j]:
            # print ('j:',j)
            # print ('key:',key)
            # if non_sorted_list[j]>key :
            non_sorted_list[j+1]=non_sorted_list[j]
            j-=1
        non_sorted_list[j+1]=key
            
    print (non_sorted_list)
    return(non_sorted_list)
    
print('*******************************************')    
my_test_list=[7,1,5,2,6,3,1,9,11,7,23]
my_sorted_list=Insertion_Sort(my_test_list)

#------------------------------------------------------------
#second function: 
#analize of web services performance in infrastructure 
#architecture (the forth layer of an enterprise architecture)
#------------------------------------------------------------

def Web_Service_Analyze(services,resp_times,exe_CPU_costs,exe_Mem_costs,exe_Disk_costs,exe_Net_costs):
    '''
    

    Parameters
    ----------
    services : list pof str
        the list of services..
    resp_times : list of int
        the list of the responce times of specified web services, measured in millisecond.
    exe_CPU_costs : list of int
        list of the percent of execution costs in case 
        of cpu for specified web services.
    exe_Mem_costs : list of int
        list of the percent of execution costs in case 
        of memory for specified web services.
    exe_Disk_costs : list of int
        list of the percent of execution costs in case
        of disk for specified web services.
    exe_Net_costs : list of int
        list of the percent of execution costs in case
        of network for specified web services.

   Returns
   a list of services with their useful information for easily analyze.

'''
    web_services_analyze_data=[]
    for i in range(len(services)):
        serv_name=services[i]
        resp_time=resp_times[i]
        exe_CPU_cost=exe_CPU_costs[i]
        exe_Mem_cost=exe_Mem_costs[i]
        exe_Disk_cost=exe_Disk_costs[i]
        exe_Net_cost=exe_Net_costs[i]
        
        
        
        if resp_time < 100:
            resp_time_st='Excellent'
        elif resp_time >= 100 and resp_time <= 200:
            resp_time_st='good'
        elif resp_time > 200 and resp_time <= 1000:
            resp_time_st='acceptable'
        elif resp_time > 1000 :
            resp_time_st='very slow'
            
            
        if exe_CPU_cost >= 0 and exe_CPU_cost < 30:
                exe_CPU_cost_st='Low Usage'
        elif exe_CPU_cost >= 30 and exe_CPU_cost <= 70:
                exe_CPU_cost_st='Optimal'
        elif exe_CPU_cost > 70 and exe_CPU_cost <= 85:
               exe_CPU_cost_st='Moderate Pressure'
        elif exe_CPU_cost > 85 and exe_CPU_cost <= 100:
                exe_CPU_cost_st='Critical'
                
        
        if exe_Mem_cost >= 0 and exe_Mem_cost < 50:
                exe_Mem_cost_st='Low Usage'
        elif exe_Mem_cost >= 50 and exe_Mem_cost <= 75:
                exe_Mem_cost_st='Optimal'
        elif exe_Mem_cost > 75 and exe_Mem_cost <= 90:
               exe_Mem_cost_st='Moderate Pressure'
        elif exe_Mem_cost > 90 and exe_Mem_cost <= 100:
                exe_Mem_cost_st='Critical'
                
                
        if exe_Disk_cost >= 0 and exe_Disk_cost < 50:
                exe_Disk_cost_st='Low Usage'
        elif exe_Disk_cost >= 50 and exe_Disk_cost <= 75:
                exe_Disk_cost_st='Optimal'
        elif exe_Disk_cost > 75 and exe_Disk_cost <= 90:
               exe_Disk_cost_st='Moderate Pressure'
        elif exe_Disk_cost > 90 and exe_Disk_cost <= 100: 
                exe_Disk_cost_st='Critical'
                
                
        if exe_Net_cost >= 0 and exe_Net_cost < 50:
                exe_Net_cost_st='Low Usage'
        elif exe_Net_cost >= 50 and exe_Net_cost <= 75:
                exe_Net_cost_st='Optimal'
        elif exe_Net_cost > 75 and exe_Net_cost <= 90:
               exe_Net_cost_st='Moderate Pressure'
        elif exe_Net_cost > 90 and exe_Net_cost <= 100: 
                exe_Net_cost_st='Critical'
                
                
       
        web_services_analyze_data.append(('Service Name: '+ serv_name,'Response Time: '+resp_time_st,'CPU Cost: '+ exe_CPU_cost_st,'Memory Cost: '+ exe_Mem_cost_st,'Disk Cost: '+exe_Disk_cost_st,'Network Cost: '+exe_Net_cost_st))
        
    print(web_services_analyze_data)
    return (web_services_analyze_data)

print('*******************************************')

my_test_resp_times=[84,60,10,57,86]
my_test_cpu_costs=[61,75,91,43,92]
my_test_mem_costs=[35,72,61,46,88]
my_test_disk_costs=[75,39,24,15,36]
my_test_net_costs=[84,61,35,75,91]
my_test_services=['LLaMA3','IMDB','freepik','RSS','translate']
my_test_analyze_data=[]
my_test_analyze_data=Web_Service_Analyze(my_test_services, my_test_resp_times, my_test_cpu_costs,my_test_mem_costs,my_test_disk_costs,my_test_net_costs)




#part 3
#convertions

#---------------------------------------------
#calculation the monthly loss of an enterprise
#---------------------------------------------

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




print('*******************************************')
test=Convert_Persion_Hours_To_Persion_Days(767,6)
print (test)

#************************ 5 Constant Numbers********************************

#Avogadro Number (1/mol)
NA=6.22*10**23

#Faraday Constant (C/mol)
F= 96485.33

#Atomic Mass Comstant
amu=1.660538*10**(-27)

#Gas Constant (J/(mol.K))
R=8.3144

#Planck Constant (J.s)
h=6.626*10**(-34)

#Molar Volume of an Ideal Gas at STP  (L/mol)
MV=22.414

#Electron Mass (kg)
me=9.109*10**(-31)

#Proton Mass (kg)
mp=1.673*10**(-27)

#Neutron Mass (kg)
mn=1.675*10**(-27)


#=============================================================================
#=============================================================================
#=============================================================================
#=============================================================================


#************************ 2 functions in chemistry*****************************


#Defect density(Delta) of crystalline structures calculates from Scherrer Equation(D) and 1/(D**2)


import math

def Defect_Density(Beta,Theta,K=0.9,Landa=1.5406):
    '''
    
    Parameters
    ----------
    Beta : Float
        Full width at half maximum (FWHM).
    Theta : Float
        Bragg's Diffraction Angle.
    K : Float, optional
        Scherrer Constant. The default = 0.9.
    Landa : Float, optional
        X-ray Wavelength. The default = 1.5406.
    

    Returns
    -------
    D : Float
        Crystallite size (nm)
    Delta : Float
        Density of defect for crystallite structures from XED crystallography calculated from reverse the Scherrer Equation.

    '''
    Theta=math.radians(Theta)
    D=(K*Landa)/(Beta*math.cos(Theta))
    Delta=1/(D**2)
    return D,Delta

Defect_Density(0.6, 90)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Absorption Coefficient (D) in electrochemical reactions calculates from Randles-Sevcik aquation in Cyclic voltametry technique.
# D calculation depends on type of electrochemical reaction (reversibility and irreversibility)

def Diffusion_Coefficient_Calculator(Peak_Current,A,C,Scan_Rate,n,is_reversible):
    '''
    

    Parameters
    ----------
    Peak_Current : Float
        Peak of current (A).
    A : Float
        Electrode area (cm**2).
    C : Float
        Concentration of electroactive species (mol/L).
    Scan_Rate : Float
    (mV/s)
    n : int
    number of transfered electron
    is_reversible: bool
        Type of reaction. 'True' is reversible, 'False' is irreversible

    Returns
    -------
    Diffusion Coefficient calculates from Randles-Sevcik aquation. type:float

    '''
    
    if is_reversible.upper()=='TRUE':
        Diffusion_Coefficient=(Peak_Current/((2.65*10**5)*(n**(3/2))*A*C*(Scan_Rate**(1/2))))**2
        
    elif is_reversible.upper()=='FALSE':
        Diffusion_Coefficient=(Peak_Current/(0.446*(n**(3/2))*A*C*(Scan_Rate**(1/2))))**2
        
    else:
        raise TypeError ('Please enter a valid answer (yes or not)') 
    
    return(Diffusion_Coefficient)
    
    
    #=============================================================================
    #=============================================================================
    #=============================================================================
    #=============================================================================
    

    #*************************** 2 Converters************************************
    
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


import math
#-------------- Part 1 ----------------------------

eps=2.2204e-16           # epsilon
e=2.71828                # adad Neper
pi=3.14
phi=1.618033            # Golden Ratio
G=6.67e-11              # Gravitational constant
g=9.82                  #shetab Geranesh

#---------------------------------------------------
#----------------- Part 2 --------------------------

#+++++++++++++++++ first +++++++++++++++++++++++++++
def Component(a,b):
    '''
    تصویر اسکالر b بر a  
(اندازه بردار تصویر)       

    Parameters
    ----------
    a : list
        multi-dimensional vector
    b : list
        multi-dimensional vector

    Returns
    -------
   c : float
       The component of b on a

    '''
    ab=0
    aa=0
    for i in range(0,len(a)):
           ab=ab+a[i]*b[i]
           aa=aa+a[i]*a[i]
    aas=math.sqrt(aa)
    c=(ab/aas)
    return c

#+++++++++++++++++++++ second +++++++++++++++++++++++++++++++
def Distance(x,y,z,a,b,c,d):
    '''
    The formula of the distance between a point and a plane in
    three-dimensional space is as follows:
        D=(abs(a*x+b*y+c*z+d))/sqrt(a**2+b**2+c**2)

    Parameters
    ----------
    x : int
        The x-axis component.
    y : int
        The y-axis component.
    z : int
        The z-axis component.
    a : int
        The coefficient of x in the plane equation.
    b : int
        The coefficient of y in the plane equation.
    c : int
        The coefficient of z in the plane equation.
    d : int
        The constant in the plane equation.

    Returns
    -------
    D : float
        The distance of point-to-plane
        
   Example : Distnace between (1,-2,4) and the plane 13*x-6*y-5*z+3=0
             is    0.527 

    '''
    if a==b==c==0:
        return None
    else:
        tt=abs(a*x+b*y+c*z+d)
        ss=math.sqrt(a**2+b**2+c**2)
        D=tt/ss
        return D

#----------------- Part 3  (Conventor) --------------------------------
def Mesghal_to_Gram(mesghal):
    '''
    Tabe tabdile mesghal be gram

    Parameters
    ----------
    mesghal : number
        DESCRIPTION.

    Returns
    -------
    gram : number
        DESCRIPTION.

    '''
    gram=mesghal*4.608
    return gram



def Gram_To_Mesghal(gram):
    mesghal=gram/4.608
    return mesghal


#--------------- This function get input and Give output
'''
    Function:
        Print Current time and passes it to main
    Input:
        boolean
    output:
        datatime    
        
'''
def Time_Func_Print(booolean):
    x= datetime.datetime.now()
    print(f'today is : {x} ')
    return x
'''
    Function:
        Get the Operation from User
    Input:
        boolean 
    output:
        int: type of operation     
        
'''
def Operation_Choosing(booolean):
    print("\nWhich sort of Operation would you like?")
    while(booolean==True):
        i=input("\n\nfor Choosing Today's past Hours press 1 \nfor Choosing Today's past Minutes press 2\nfor Choosing Today's past second press 3 \nfor Choosing the amount of days has passed since today press 4 : ")
        if( int(i)==1 or int(i)==2 or int(i)==3 or int(i)==4):
            booolean=False
            break
        print("You have chosen wrong operation! please try again!")
            
    return int(i)
    

'''
    Function:
        Convertion 
    Input:
        boolean 
    output:
        int hour
        int minute
        int second
        int type of operation     
'''
def Unit_conversion_High_to_Low(const,hour,minute,sec,operation):   
    if(operation==2):
       return hour*const+minute
    if(operation==3): 
        temp = (hour*const+((minute*const)/60))+int(sec)
        return temp
    else: 
        return -1
'''
    Function:
        printing result
    
    Input:
        int operation
        int hour
        int minute
        int secound
        int converted time
    output:
        void     
'''    
def Print_result(operation,hour,minute,second,temp):   
    if(operation==1):
        if(hour<=1):
            print(f'\nCurrent time is: {hour}:{minute}:{second}')
            print(f'\nUp to now {hour} hour has passed')
        elif(hour>1):
            print(f'\nCurrent time is: {hour}:{minute}:{second}')
            print(f'\nUp to now {hour} hours has passed')
    elif(operation==2):
        if(minute<=1):
            print(f'\nCurrent time is: {hour}:{minute}:{second}')
            print(f'\nUp to now {temp} minute has passed')
        elif(minute>1):
            print(f'\nCurrent time is: {hour}:{minute}:{second}')
            print(f'\nUp to now {temp} minutes has passed')
    elif(operation==3):
        if(second<=1):
            print(f'\nCurrent time is: {hour}:{minute}:{second}')
            print(f'\nUp to now {temp} second has passed')
        elif(second>1):
            print(f'\nCurrent time is: {hour}:{minute}:{second}')
            print(f'\nUp to now {temp} seconds has passed')
    elif(operation==4):
        if(temp<=1):
           print(f'\nUp to now {temp} day has passed')
        elif(hour>1):
           print(f'\nUp to now {temp} days has passed')
    
    
'''
    Function:
        convertion: Counting Days
    
    Input:
        
        int mounth
        int day
    output:
        int passed day    
'''    
def Unit_Convertion_Day(mounth,day):
    arr=[['jan',31],['feb',28],['mar',31],['apr',30],['may',31],['june',30],['july',31],['agust',31],['september',30],['october',31],['november',30],['december',31]]
    temp=0
    for i in range(mounth-1):
        temp+=arr[i][1]
    return temp+day-1

#--------------- Constant Numbers  
hour_to_minute=60
minute_to_sec=60
hour_to_sec=3600
gama_const=0
b=True 
current_time = Time_Func_Print(b)
operation=Operation_Choosing(b)
if(operation==1):
    Print_result(1,current_time.hour, current_time.minute, current_time.second,0)
elif(operation==2): #converting hour to minute
    temp = Unit_conversion_High_to_Low(hour_to_minute,current_time.hour,current_time.minute,current_time.second,operation)
    Print_result(2,current_time.hour, current_time.minute, current_time.second,temp)
elif(operation==3): #converting hour to second
    temp = Unit_conversion_High_to_Low(hour_to_sec,current_time.hour,current_time.minute,current_time.second,operation)
    Print_result(3,current_time.hour,current_time.minute,current_time.second,temp)
elif(operation==4):
    temp=Unit_Convertion_Day(current_time.month,current_time.day)
    Print_result(4,current_time.hour,current_time.minute,current_time.second,temp)


#----------1----------
N=6.02214076*(10**23) #Avagadro's Number
R=8.3145              #Gas Constant
K=1.380649*(10**-23)  #Boltzmann constant
e=2.71828             #Euler's Number

Z_A=0.01              #Flexure Testing of Polymers- A:rate of straining of outer surface at 5% strain 
Z_b=0.1               #Flexure Testing of Polymers- B:rate of straining of the outer surface intended for materials that may not break at 5% strain


#----------2----------
def Change_In_Pressure_Of_a_Mercury_OR_Water_Column_Manometer(Type_Of_Fluid,h,g=9.81):
    
 if Type_Of_Fluid.upper()=='water':
    Delta_P=1*h*g
    return Delta_P
     
 if Type_Of_Fluid.upper()=='mercury':
    Delta_P=13.6*h*g
    return Delta_P   

 else:
     raise ValueError ("Invalid Type of Fluid")


def Mass_Of_Rubber_In_Internal_Mixer(filler_percentage,type_of_rotors,V):
    if filler_percentage>65:
           if type_of_rotors=='tangential':
            M=1.3*0.8*V
            return M
           
           if type_of_rotors=='intermix':
            M=1.3*0.6*V
            return M
        
    if filler_percentage<=65:
           if type_of_rotors=='tangential':
            M=1*0.8*V
            return M
           
           if type_of_rotors=='intermix':
            M=1*0.6*V
            return M
        
        
#----------3----------
def Kg_to_Lbm(Kg):
    Lbm=0.4535*Kg
    return Lbm
   

def Lbm_to_Kg(Lbm):
    Kg=2.20462*Lbm
    return Kg






###########################################################################
'''Section 1:
    Constant numbers in engineering'''
###########################################################################

# Ideal gas constant
R = 8.314                # J/mole.K


#Avogadro's number
NA = 6.023*10**23

# Faraday's Constant
F= 96485                #C/mole


#Boltzmann Constant
k_B = 1.381*10**-23         #J/K



#Planck's Constant
h = 6.626*10**-34           #JS




###########################################################################
'''Section 2:
    Convertor function from atmosphere to mmHg and vice versa'''
###########################################################################

#  atmosphere to mmHg 

def Atmosphere_MmHg (P):
    '''
    Parameters
    ----------
    P : float
        pressure in atm unit.

    Returns
    pressure from atm to mmHg
    '''
   
    P_mmHg= P*760
    return P_mmHg
    
# mmHg to atm

def MmHg_Atmosphere (P):
    '''
    Parameters
    ----------
    P : float
        pressure in mmHg unit.

    Returns
    pressure from mmHg to atm

    '''
    P_atm= P/760
    return P_atm



###########################################################################
'''Section 3-1

    A function that takes k and n and returns the list of k prime numbers after n'''
###########################################################################
    
def Addad_Aval (n,k):
    '''
    Parameters
    ----------
    n : int
        Desired number.
    k : int
        Number of prime numbers after n.

    Returns
    -------
    list_adad : list
        lists  the "k" prime numbers after "n" .
    '''
        
    list_adad=[]
    count=0
    while count<k:
        if (n+1)%2 ==0 and (n+1)!=2:
            n=n+1
        elif (n+1)%3 ==0 and (n+1)!=3:
            n=n+1
        elif (n+1)%5 ==0 and (n+1)!=5:
                n=n+1
        elif (n+1)%7 ==0 and (n+1)!=7:
            n=n+1
        else:
            n=n+1
            count=count +1
            list_adad.append (n)
    return (list_adad)





###########################################################################
'''Section 3-2

    A function that takes the list of grades of a student and the 
    number of units of each course and return the averag of the student. 
    Lessons with a score below 10 are not included in the average and their 
   numbers printe as the number of ineffective lessons'''
###########################################################################
def Average (Number, Vahed):
      Sum=0
      count1=0
      count2=0
      for i in range (0, len(Number)):
          if Number[i]>=10:
              Sum=Sum+(Number[i]*Vahed[i])
              count1=count1+1
          else:
              count2=count2+1
          mean= Sum/count1
      print('There are', count2, 'non-effective number')
      return mean      






###########################################################################
'''Section 3-3

    A function that Calculate the output concentration and flowrate from
    a tank where two flows with different concentrations and flow rates
    enter it'''
###########################################################################


def Concentration_Calculator (CA1, unit_CA1, CA2, unit_CA2, Q1, unit_Q1, Q2, unit_Q2):
    '''
    Parameters
    ----------
    CA1 : float
        Concentratoion of the first flow
    unit_CA1: str
        Unit of concentratoion of the first flow
    CA2 : float
        DESCRIPTION.Concentratoion of the second flow
    unit_CA2: str
        Unit of concentratoion of the second flow
    Q1 : float
        Molar volume of the first flow
    unit_Q1: str
        Unit of molar volume of the first flow
    Q2 : float
        Molar volume of the second flow
    unit_Q2: str
        Unit of molar volume of the second flow

    Returns
    -------
    Flowrate and concentration of the output flow
    '''
    
    if unit_CA1== 'mole/lit':
        CA1= CA1*1000
    elif unit_CA1== 'mole/ml':
        CA1= CA1*1000000
        
    if unit_CA2== 'mole/lit':
        CA2= CA2*1000
    elif unit_CA2== 'mole/ml':
        CA2= CA2*1000000
        
    if unit_Q1== 'lit':
        Q1= Q1/1000
    elif unit_Q1== 'ml' or unit_Q1== 'cc':
        Q1= Q1/1000000 
        
    if unit_Q2== 'lit':
            Q2= Q2/1000
    elif unit_Q2== 'ml' or unit_Q2== 'cc':
            Q2= Q2/1000000 
            
    V_out = Q1+Q2
    C_out= (CA1*Q1 + CA2*Q2)/V_out
    return V_out
    return C_out




###########################################################################
'''Section 3-4

    A function that Determines the solubility of substances at different 
    temperatures using  of the constants related to that substance and 
    the end of the calculations when the concentration is greater than 1'''
###########################################################################


def Solubility (a,b,c,T, T_unit, T_end):
    '''
     Parameters
    ----------
    a : float
        Constant.
    b : float
        Constant.
    c : float
        Constant.
    T : float
        Temperature.
    T_unit:str
        Unit of temperature.
    T_end: float
        Final temperature

    Returns
    -------
    Determination solubility at different temperatures upto solubility< 1
    '''
#تبدیل دما به کلوین

    if T_unit=='C':
        T_unit=T_unit+273.15
    elif T_unit=='F':
        T_unit= ((T_unit-32)/1.8)+273.15
    elif T_unit=='R':
        T_unit=T_unit*5/9
        
# محاسبه حلالیت و توقف محاسبات با رسیدن به غلظت بالاتر از 1  
  
    for T in range (T,T_end):
        S=a*T**2+b*T+c
        if S<1:
            continue
        else:
            break
    return S
        
'''-----------------------------------------------------------------------
------------------------------------------------------------------------'''   

#============================================================================#
#============================================================================#
#===========================A1_Mohammad_Shafiee==============================#
#============================================================================#
#============================================================================#

#---------------------------------<Part 1>-------------------------------------
Electron_Charge_Constant=1.602*(10**-19)
Universal_Gas_Constant=8.314
Avogadro_Num=6.022*(10**23)
Planck_Constant=6.62607*(10**-34)
Young_Modulus_Steel=200

#==============================================================================
#---------------------------------<Part 2>-------------------------------------


#-------->Function 1: Tresca Criterion<---------

#-------->In the case that the principal stresses are given
def Tresca_Yield_For_Principal_Stresses(c,hardness,sigma_1,sigma_2,sigma_3,/):
    '''
    

     Parameters
     ----------
     c : float
         C is a coefficient that is between 0.3 and 0.4 depending on the type of material.
     hardness : float
         brinell and vickers hardness.
     sigma_y : float
         Uniaxial tensile yield strength of material.(Mpa)
     sigma_1 : float
         Principal Stresses.
     sigma_2 : float
         Principal Stresses.
     sigma_3 : float
         Principal Stresses.

     Returns
     -------
     Tresca : float
         Tresca Yield(Mpa)

    '''
    sigma_y = c*hardness 
    
    Max=max(sigma_1,sigma_2,sigma_3)
    Min=min(sigma_1,sigma_2,sigma_3)
    
    Tresca=Max-Min
    
    if Tresca < sigma_y :
        print("It is in the elastic area")
    elif Tresca==sigma_y :
        print("It is on the verge of plastic deformation")
    
    elif Tresca > sigma_y :
        print("Plastic deformation has occurred")
        
    return Tresca
        
#==============================================================================   
#--------->If we have biaxial tension
def Tresca_Yield_For_Biaxial_Tension(c,hardness,sigma_xx,sigma_yy,tau_xy,/):
    '''
    

    Parameters
    ----------
    c : float
        C is a coefficient that is between 0.3 and 0.4 depending on the type of material.
    hardness : float
        brinell and vickers hardness.
    sigma_y : float
        Uniaxial tensile yield strength of material.(Mpa)
    sigma_y : float
        Uniaxial tensile yield strength of material
    sigma_xx : float
        
    sigma_yy : float
        
    tau_xy : float
        

    Returns
    -------
    Tresca : TYPE
        DESCRIPTION.

    '''
    import math
    sigma_y = c*hardness
    sigma_max=(sigma_xx + sigma_yy)/2 + math.sqrt(((sigma_xx-sigma_yy)/2)**2 + tau_xy**2)
    sigma_min=(sigma_xx + sigma_yy)/2 - math.sqrt(((sigma_xx-sigma_yy)/2)**2 + tau_xy**2)

    Tresca=sigma_max-sigma_min
    
    if Tresca < sigma_y :
        print("It is in the elastic area")
    elif Tresca==sigma_y :
        print("It is on the verge of plastic deformation")
    elif Tresca > sigma_y :
        print("Deformation is plastic")
            
    return Tresca
    
#==============================================================================
#-------->Function 2: Total Solidification Time of Casting


def Total_Solidification_Time_of_Casting(volume,surface_area,Cm,n,/):
    '''
    

    Parameters
    ----------
    volume : float
        DESCRIPTION.
    surface_area : float
        DESCRIPTION.
    Cm : float
        It's amount varies depending on the type of material
    n : float
        It's amount varies depending on the type of material but is 
        usually equal to 2
        
    Returns
    -------
    Total_Solidification_Time : float

    '''

    
    Total_Solidification_Time = Cm*((volume/surface_area)**2)
    return Total_Solidification_Time

#==============================================================================
#--->Function 3: Miner's Rule for Cumulative Damage Theory and Fatigue Failures
def Miner_Rule_Fatigue(n,N):
    '''
    

    Parameters
    ----------
    n : list
        number of cycles at stress level 

    N : list
        number of cycles to failure at stress level

    Returns
    -------
    sigma : float

    '''
    sigma=0
    for i in range(0,len(n)):
        for j in range(0,len(N)):
            if i==j:
                sigma=sigma+(n[i]/N[j])
    
    if sigma<1:
        print("There is no fatigue yet")
    elif sigma>=1:
        print("Failure occurs at least according to Miner's rule")
                    
    return sigma

#==============================================================================
#---------------------------------<Part 3>-------------------------------------
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



#Part 1 ______________________________________________________________________________________________

E=190          #Young's modulus of steel at room temperature (GPA).
e=1            #Coefficient of restitution for perfectly elastic collision.
K=2.1          #Bulk's modulus of water at room temperature (GPA).
air_p=1013.25  #Air pressure at sea level (hPA).
g=9.81         #The acceleration due to gravity, Near Earth's surface (m/s2).

#Part 2_______________________________________________________________________________________________

import math

def Indeterminate_degree_of_truss(m,j):
    ''' 
    This function calculates the degree of indeterminacy of the truss by taking 'm' number of truss members and 'j' number of nodes.
    
    Parameters
    ----------
    m : int
        Number of members truss.
    j : int
        Number truss node members.

    Returns
    -------
    n : float
        The indeterminate degree of truss.
        
        if n=0 : The truss is determinate.
        if n>0 : The truss is statically unstable.
        if n<0 : The truss is indeterminate.
        
    '''
    Equations=2*j
    Unknowns=3+m
    global n
    if (Equations==Unknowns):
        n=0
    if Equations>Unknowns:
        n=Equations-Unknowns
    if Equations<Unknowns:
        n=Equations-Unknowns
    return n
    

def Friction_Law(mass,angle,/,MOTCS,moving):
    '''
    This function calculates the friction force on the object by taking the mass of the desired object, the material of the two contact surfaces, the angle of the contact surface with the horizon and the movement state of the object.
    
    Parameters
    ----------
    mass : int
        The mass of the desired object (KG).
    angle : int
        The angle of the contact surface with respect to the horizon (RAD).
    MOTCS : str
        Material of two contact surfaces.
        The contact surfaces considered :
        STEEL ON STEEL
        STEEL ON ALUMINUM
        STEEL ON COPPER
        COPPER ON CAST IRON
        COPPER ON GLASS
        GLASS ON GLASS
        RUBBER ON DRY CONCRETE
        RUBBER ON WET CONCRETE
        TEFLON ON TEFLON
    moving : str
        Is the desired object moving? (YES OR NO).

    Returns
    -------
    F : float
        Frictional force applied to the target object (N).

    '''
    global F
    m=moving.upper()
    M=MOTCS.upper()
    if m=='YES':     
        if M=='STEEL ON STEEL':
            F=mass*9.81*math.cos(angle)*0.57
        if M=='STEEL ON ALUMINUM':
            F=mass*9.81*math.cos(angle)*0.47
        if M=='STEEL ON COPPER': 
            F=mass*9.81*math.cos(angle)*0.36
        if M=='COPPER ON CAST IRON':
            F=mass*9.81*math.cos(angle)*0.29
        if M=='COPPER ON GLASS':
            F=mass*9.81*math.cos(angle)*0.53
        if M=='GLASS ON GLASS':
            F=mass*9.81*math.cos(angle)*0.40
        if M=='RUBBER ON DRY CONCRETE':  
            F=mass*9.81*math.cos(angle)*0.8
        if M=='RUBBER ON WET CONCRETE':   
            F=mass*9.81*math.cos(angle)*0.25
        if M=='TEFLON ON TEFLON':   
            F=mass*9.81*math.cos(angle)*0.04
        return F
    if m=='NO':
        if M=='STEEL ON STEEL':
            F=mass*9.81*math.cos(angle)*0.74
        if M=='STEEL ON ALUMINUM':
            F=mass*9.81*math.cos(angle)*0.61
        if M=='STEEL ON COPPER': 
            F=mass*9.81*math.cos(angle)*0.53
        if M=='COPPER ON CAST IRON':
            F=mass*9.81*math.cos(angle)*1.05
        if M=='COPPER ON GLASS':
           F=mass*9.81*math.cos(angle)*0.68
        if M=='GLASS ON GLASS':
           F=mass*9.81*math.cos(angle)*0.94
        if M=='RUBBER ON DRY CONCRETE':  
           F=mass*9.81*math.cos(angle)*1
        if M=='RUBBER ON WET CONCRETE':   
           F=mass*9.81*math.cos(angle)*0.3
        if M=='TEFLON ON TEFLON':   
           F=mass*9.81*math.cos(angle)*0.04
        return F
        
#Part 3_______________________________________________________________________________________________

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


#=================================
#================================
#===========================
#====================
#=============
#=======
#=======
#============
#====================
#=========================
#=================================


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



Data=pd.read_excell('C:\\Users\amina\Desktop\LL\signal-to-noise ratio')



def Signal_To_Noise_Ratio(data,application):
    '''

    Parameters
    ----------
    data : DataFrame
        consists of your experimental data in 3 columns:
            1- location: the place you've measured signal and noise
            2- signal strength: the power of signal in dbm
            3- noise power: the power of noise in dbm
    application : string
        there is 3 application available:
            1- plot signal: plots signal column
            2- plot noise: plots noise column
            3- plot SNR: plots the signal-to-noise ratio

    Returns
    -------
    mx : float
        the maximum signal-to-noise ratio in dbm

    '''
    location=np.array(data['location'])
    signal=np.array(data['Signal Strength'])
    noise=np.array(data['Noise Power'])
    snr=signal-noise
    
    
    if str(application).lower()=='plot signal' :
        plt.plot(location,signal)
        plt.title('signal power in every place')
        plt.xlabel('place')
        plt.ylabel('signal power in dbm')
        plt.grid()
        plt.show()       
    elif str(application).lower()=='plot noise' :
        plt.plot(location,noise)
        plt.title('noise power in every place')
        plt.xlabel('place')
        plt.ylabel('noise power in dbm')
        plt.grid()
        plt.show()
    elif str(application).lower()=='plot snr' :
        plt.plot(location,snr)
        plt.title('SNR in every place')
        plt.xlabel('place')
        plt.ylabel('SNR in db')
        plt.grid()
        plt.show()
        
        
    mx=snr.max()
    return mx




#climate_change


import numpy as np

import matplotlib.pyplot as plt

import pandas as pd

a=np.random.uniform((54<62.276),(15.34<16.82),size=(6,3))
data=pd.DataFrame(a,index=['year1980','year1990','year2000','year2010','year2020','year2024'],columns=['Fahrenheit','Celsius','Average Temperature'])
data.max()

data.dropna(inplace=True)
data.info()

plt.hist(a)
plt.xlabel('Celsius')
plt.ylabel('Fahrenheit')
plt.show()
x=['1980','1990','2000','2010','2020','2024']
y=np.array([0,10,20,25,30,35,40])

#Fahrenheit_Celsius

def Fahrenheit_Celsius(a):
    Fahrenheit=np.array(a['Celsius'])
    Celsius=np.array(a['Fahrenheit'])
    
    if a=='plot':
        plt.plot(Fahrenheit,Celsius)
        plt.title('Fahrenheit','Celsius')
        plt.show()
        
    elif a=='Average Temperature in 1980':
        Average_Temperature1980=(59.62 + 15.344444)/2
        return Average_Temperature1980
    
    elif a=='Average Temperature in 1990':
        Average_Temperature1990=(59.8 + 15.44444)/2
        return Average_Temperature1990
    
    elif a=='Average Temperature in 2000':
        Average_Temperature2000=(54 + 12.2222)/2
        return Average_Temperature2000
    
    elif a=='Average Temperature in 2010':
        Average_Temperature2010=(58.4 + 14.66667)/2
        return Average_Temperature2010
    
    elif a=='Average Temperature in 2020':
        Average_Temperature2020=(58.82 + 14.9)/2
        return Average_Temperature2020
    
    elif a=='Average Temperature in 2024':
        Average_Temperature2024=(62.276 + 16.82)/2
        return Average_Temperature2024

def climate_change(data,application):
    '''
    

    Parameters
    ----------
    data : climate change temperature
        <class 'pandas.core.frame.DataFrame'>
        Index: 6 entries, year1980 to year2024
        Data columns (total 3 columns):
         #   Column               Non-Null Count  Dtype  
        ---  ------               --------------  -----  
         0   Fahrenheit           6 non-null      float64
         1   Celsius              6 non-null      float64
         2   Average Temperature  6 non-null      float64
        dtypes: float64(3)
        memory usage: 192.0+ bytes
    application : Fahrenheit to Celsius
        The formula for converting Fahrenheit to Celsius is C = 5/9(F-32).
        Fahrenheit and Celsius are the same at -40°. At ordinary temperatures, Fahrenheit is a larger number than Celsius.

    climate_change:
    the combined land and ocean temperature has increased at an average rate of 0.11° Fahrenheit (0.06° Celsius) per decade since 1850,
    or about 2° F in total. The rate of warming since 1982 is more than three times as fast: 0.36° F (0.20° C) per decade.
    '''




import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

#%%
def polarization_control(Data,Application):
    '''
    

    Parameters
    ----------
    Data : DataFrame : The results of the polymerization of a polymer control laboratory are given 
                       that included four columns: time, temperature, pressure, 
                       and reaction percentage of the polymer product.
        DESCRIPTION.
    Application : float
        Six applications are done in this function.
        1) Application = 'temp_time' : Temperature is plotted according to time. 
        2) Application = 'pressure_time' : Pressure is plotted according to time. 
        3) Application = 'Percent_time' : The percentage of reaction is plotted over time.
        4) Application = '100% reaction': If the percentage of polymerization reaction proceed to 100, the temperature and pressure of polymerization is printed and returned.
        5) Application = 'Max_pressyre': It returns maximum pressure of process.
        6) Application = 'Max_temp': It returns maximum temperature of process.
    Returns
    -------
    float
        It returns temperature and pressure of process according to related application.

    '''
    
    time = np.array(Data['time'])
    temp = np.array(Data['temp'])
    pressure = np.array(Data['pessure'])
    reaction_percent = np.array(Data['percent'])
    
    if Application == 'temp_time':
        plt.plot(time, temp, c = 'g',linewidth = 1.5)
        xylable_font={'family': 'serif',
                'color': 'black' ,
                'size': 16 }
        title_font={'family': 'serif',
                'color': 'black' ,
                'size': 16 }
        plt.title('Temperature variation',fontdict = title_font)
        plt.xlabel('time(s)',fontdict = xylable_font)
        plt.ylabel('Temperature (C)',fontdict = xylable_font)
        # plt.legend(['Temperature'])
        plt.show()

    elif Application == 'pressure_time':
        plt.plot(time, pressure , c = 'r',linewidth = 1.5)
        xylable_font={'family': 'serif',
                'color': 'black' ,
                'size': 16 }
        title_font={'family': 'serif',
                'color': 'black' ,
                'size': 16 }
        plt.title('Pressure variation',fontdict = title_font)
        plt.xlabel('time(s)',fontdict = xylable_font)
        plt.ylabel('Pressure (Pa)',fontdict = xylable_font)
        # plt.legend(['Pressure'])
        plt.show()

    elif Application == 'Percent_time':
        plt.plot(time, reaction_percent, c = 'b',linewidth = 1.5)
        xylable_font={'family': 'serif',
                'color': 'black' ,
                'size': 16 }
        title_font={'family': 'serif',
                'color': 'black' ,
                'size': 16 }
        plt.title('Progress variation',fontdict = title_font)
        plt.xlabel('time(s)',fontdict = xylable_font)
        plt.ylabel('Reaction percent',fontdict = xylable_font)
        # plt.legend([ 'reaction_percent'])
        plt.show()

    elif Application == '100% reaction':    
        reaction_percent_p=np.arange(0,301)/3
        L = len(reaction_percent_p)
        for i in range(L):
            if reaction_percent_p[i] == 100:
                print('tempreature and pessure for 100% progress are','tempreature =',temp[i],'(C)', 'pessure=', pressure[i],'(Pa)')
                return   (temp[i], pressure[i])     
    elif Application == 'Max_pressyre':    
         return pressure.max()
    elif Application == 'Max_temp':    
         return temp.max()
#%% Calling function
time = np.reshape(np.arange(0,301),(301,1))
temperature= np.reshape(np.random.normal(loc=40,scale=1,size=(301)),(301,1))
pressure = np.reshape(np.random.normal(loc=2,scale=0.2,size=(301)),(301,1))
reaction_percent= np.reshape(np.arange(0,301)/3,(301,1))

data = np.concatenate((time,temperature,pressure,reaction_percent),axis=1)
Data = pd.DataFrame(data,columns=['time','temp','pessure','percent'])


polarization_control(Data,'temp_time')


Data=pd.read_excel('/Users/elnazmac/Desktop/Desulfurization.xlsx')

#################FUNCTION########################################

def Desulfurization_Rate(Data,application):
    '''

    Parameters
    ----------
    Data : Data Frame
        experimental data (excel).
    application : 
        1.plot
        2.Max_Removal_With_Ultrasonic
        3.Max_Removal_Without_Ultrasonic

    Returns

    '''
    x=np.array(Data['Time'])
    y1=np.array(Data['Desulfurization_With_Ultrasonic'])
    y2=np.array(Data['Desulfurization_Without_Ultrasonic'])
    
    if application=='plot':
        plt.plot(x,y1,marker='*',mec='r',mfc='y',ms=14,ls='-',linewidth=5,color='r',label='With_Ultrasonic')
        plt.plot(x,y2,marker='*',mec='g',mfc='y',ms=14,ls='--',linewidth=5,color='g',label='Without_Ultrasonic')
        
        myfont={ 'family': 'serif'   ,
                'color':  'red'  ,
                'size':  15   }
        
        plt.xlabel('Time')
        plt.ylabel('Desulfurization')
        plt.title('Sulfur_Removal_Plot',fontdict=myfont)
        plt.legend()
        plt.show()
    
    elif application=='Max_Removal_With_Ultrasonic':
        Max_Removal_With_Ultrasonic=y1.max()
        return Max_Removal_With_Ultrasonic
    
    elif application=='Max_Removal_Without_Ultrasonic':
        Max_Removal_Without_Ultrasonic=y2.max()
        return Max_Removal_Without_Ultrasonic



#karbord tabe1
Data=pd.read_excel('/Users/elnazmac/Desktop/Desulfurization.xlsx')
Desulfurization_Rate(Data,'plot')

#karbord tabe2
Data=pd.read_excel('/Users/elnazmac/Desktop/Desulfurization.xlsx')
Desulfurization_Rate(Data,'Max_Removal_With_Ultrasonic')  

#karbord tabe3
Data=pd.read_excel('/Users/elnazmac/Desktop/Desulfurization.xlsx')
Desulfurization_Rate(Data,'Max_Removal_Without_Ultrasonic')     
    


====================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
XRD=pd.read_excel('C:/Users/Nikan/Desktop/ZnO.xlsx')
def XRD_ZnO(XRD,application):
    '''
    

    Parameters
    ----------
    XRD : DataFrame
        Data containing XRD data.
    application : str
        Type of application 'plot','FWHM','Scherrer'.
        plot:To draw the figure.
        FWHM:Width at Half Maximum.
        Scherrer:To calculate the crystallite size.

    Returns
    FWHM,Scherrer
    -------
    None.

    '''
    Angles=np.array(XRD['Angle'])
    Intensities=np.array(XRD['Det1Disc1'])
    if  application=='plot':
        plt.plot(Angles,Intensities,c='red')
        plt.title('XRD Pattern')
        plt.xlabel('2theta (degrees)')
        plt.ylabel('Intensity')
        plt.show()
    elif application in ['FWHM', 'Scherrer']:
        max_intensity = np.max(Intensities)
        half_max = max_intensity / 2
        indices = []
        half_max = max_intensity / 2
        for i in range(len(Intensities)):
           if Intensities[i] >= half_max:
               indices.append(i)

        
        if len(indices) > 0:
            left_index = np.min(indices)
            right_index = np.max(indices)
       
    
            FWHM = Angles[right_index] - Angles[left_index]
            if application == 'FWHM':
                return FWHM
           
            elif application =='Scherrer':
                mean_2theta = Angles[indices].mean()


                theta = mean_2theta / 2
                FWHM_rad = ((3.14/180)*FWHM)
                theta_rad = ((3.14/180)*theta)  
                crystal_size = (0.9 * 1.5406) / (FWHM_rad * np.cos(theta_rad))
                
                return crystal_size


crystal_size = XRD_ZnO(XRD, 'Scherrer')


XRD_ZnO(XRD,'plot')
fwhm = XRD_ZnO(XRD, 'FWHM')

crystal_size = XRD_ZnO(XRD, 'Scherrer')
print(fwhm,crystal_size)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_excel(r'C:\Users\Admin\Downloads\Telegram Desktop\XRD.xlsx')

def XRD(data,application):
    '''
    This function plots the XRD curve .

    Parameters
    ----------
    data : DataFrame
        data is XRD data include Intensity (a.u.) and 2θ (degree).
    application : str
        application is the function that you want to apply to your data:
            - plot : ُPlot the XRD curve
            - maxintensity : Determining the maximum intensity.
            - meantheta : Determining the mean of theta angles.
            

    Returns
    -------
    plot, maxintensity, meantheta

    '''
    
    data.columns = ('2θ (degree)','Intensity (a.u.)')
    Intensity = np.array(data['Intensity (a.u.)'])
    Theta = np.array(data['2θ (degree)'])
    
    if application.upper() == 'PLOT':
        font1 = {'family':'Times New Roman', 'color':'black', 'size':16}
        font2 = {'family':'Times New Roman', 'color':'black', 'size':14}
        plt.plot(Theta, Intensity, linewidth=0.8, c='k')
        plt.title('XRD', fontdict = font1)
        plt.xlabel('2θ (degree)', fontdict = font2)
        plt.ylabel('Intensity (a.u.)', fontdict = font2)
        plt.show()
        
    elif application.upper() == 'MAXINTENSITY':
        maxintensity = Intensity.max()
        return maxintensity
        
    elif application.upper() == 'MEANTHETA':
        E = 0
        for i in Theta:
            theta = i / 2
            E = E + theta        #sum of numbers
        M = E / len(Theta)       #Mean
        return M
        


XRD(data,'PLOT')


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openpyxl 

def Stress_Strain(data,application):
    '''
    this function converts F and dD to Stress and Strain by thickness(1.55mm), width(3.2mm) and parallel length(35mm).

    Parameters
    ----------
    data : DataFrame
        this DataFrame contains F(N) and dD(mm) received from the tensil test machine.
    application : str
        application determines the expected output of Stress_Strain function.

    Returns
    -------
    int, float or plot
        return may be elongation at break, strength or a plot.

    '''
    
    stress=np.array([data['F']/(1.55*3.2)])
    strain=np.array([(data['dD']/35)*100])
    if application.upper()=='ELONGATION AT BREAK':
        elongation_at_break=np.max(strain)
        print(elongation_at_break,'%')
        return elongation_at_break
    elif application.upper()=='STRENGTH':
        strength=np.max(stress)
        print(strength,'N/mm2')
        return strength
    elif application.upper()=='PLOT':
        myfont_title={'family':'sans-serif',
                      'color':'black',
                      'size':20}
        myfont_lables={'family':'Tahoma',
                       'color':'green',
                       'size':16}
        plt.plot(strain,stress,ls='--',c='g',linewidth=10)
        plt.title('Stress-Strain',fontdict=myfont_title)
        plt.xlabel('Strain(%)',fontdict=myfont_lables)
        plt.ylabel('Stress(N/mm2)',fontdict=myfont_lables)
        plt.show()
#****************function test******************        
titr=pd.DataFrame([[1,20],[2,34],[3,45],[4,67],[5,70],[4,89]],columns=['F','dD'])
Stress_Strain(titr,'plot')
Stress_Strain(titr,'elongation at break')
Stress_Strain(titr,'strength')


data_base=pd.read_excel('F-dD Data.xlsx')
tensile=pd.DataFrame(np.array(data_base),columns=['F','dD'])

#**********************clearing data***********************
tensile.info()
#out put: <class 'pandas.core.frame.DataFrame'>
#RangeIndex: 3890 entries, 0 to 3889
#Data columns (total 2 columns):
 #   Column  Non-Null Count  Dtype  
#---  ------  --------------  -----  
# 0   F       0 non-null      float64
# 1   dD      0 non-null      float64
#dtypes: float64(2)
#memory usage: 60.9 KB
# no empty cell and wrong format
count1=0
for i in tensile.index:
    if tensile.loc[i,'F']<0:
        count1=count1+1
        print(count1) #out put:0 -------> no wrong data in F column

count2=0
for j in tensile.index:
    if tensile.loc[j,'dD']<0:
        count2=count2+count1
        print(count2) #out put:0 -------> no wrong data in dD column
        
Stress_Strain(tensile,'elongation at break')
Stress_Strain(tensile,'strength')
Stress_Strain(tensile,'plot')


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
def Imerssion_Test(data,application):
    '''
    

    Parameters
    ----------
    data : .excel .csv
    columns name:time, Mg, Mg_H, Mg_Pl, Mg_HPl 
    application :
        plot:drawing the changes of weight(%) on time(days)
        More_Bioactive: the sample with more weight gain in result more bioactive
        Less_Bioactive: the sample with more weight loss in result less bioactive

    '''
    x=np.array(data['time'])
    y1=np.array(data['Mg_HPl'])
    y2=np.array(data['Mg_H'])
    y3=np.array(data['Mg_Pl'])
    y4=np.array(data['Mg'])
    if application=='plot':
        plt.plot(x,y1,marker='o',label='Mg_HPl')
        plt.plot(x,y2,marker='*',label='Mg_H')
        plt.plot(x,y3,marker='^',label='Mg_Pl')
        plt.plot(x,y4,marker='+',label='Mg')
        plt.title('The graph of changes in the weight of the samples in the SBF solution',c='r')
        plt.xlabel('Imerssion Time(day)',c='g')
        plt.ylabel('Weight Gain(%)',c='g')
        plt.legend()
        plt.show()
    elif application=='More_Bioactive':
        max_weight_gain=data[['Mg','Mg_H','Mg_Pl','Mg_HPl' ]].max()
        max_weight_gainn=max_weight_gain.max()
        more_bioactive=data.columns[data.isin ([max_weight_gainn]).any()]
        return more_bioactive
    elif application=='Less_Bioactive':
        max_weight_loss=data[['Mg','Mg_H','Mg_Pl','Mg_HPl' ]].min()
        max_weight_losss=max_weight_loss.min()
        less_bioactive=data.columns[data.isin ([max_weight_losss]).any()]
        return less_bioactive

data=pd.read_excel('weight gain.xlsx') 
Imerssion_Test(data,'Less_Bioactive')      


mport numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math 


time=np.arange(0,101).reshape(-1,1)
temp=np.random.uniform(0,101,size=101).reshape(-1,1)
pressure=np.random.uniform(0,5,size=101).reshape(-1,1)
conv=np.random.uniform(0,91,size=101).reshape(-1,1)
c=np.concatenate([time,temp,pressure,conv],axis=1)
data_dasti=pd.DataFrame(c,columns=['time','temp','pressure','conv'])
sotune4=np.array(data_dasti['conv'])
sotune2=np.array(data_dasti['temp'])

#data=pd.read_excel(r"C:\Users\asus\Desktop\A2_Prjct.xlsx.xlsx")
#dat=pd.read_csv(r'C:\Users\asus\Desktop\A2_Prjct.xlsx.xlsx')

def Conversion(data,app):
    '''
    This program is related to a chemical reaction laboratory, which has been measured in a table at a certain temperature,
    pressure and time, with a one-second interval,and gives you the highest conversion percentage at a given temperature and pressure.
    It also draws graphs of conversion percentage, temperature and pressure in terms of time.

    Parameters
    ----------
    data : DataFrame of pandas or array of numpy
        please be careful about the inputs: your table should contain about 100 index and 4  columns.
    app : str
       Only write "PLOT_TEMP" if you want to draw the figure tempreturre over time,else if you want to draw pressure on time write
       'PLOT_pressure' or write 'PLOT_CONVERSION' if you want the conversion on time figure.
       or write "MAXIMUM CONVERSION" if you want the maximum number of conversions at a specific temperature and pressure.
       Otherwise, you will see the error message below.
   
    TypeError
       The datas or application is not entered correctly

    Returns
    -------
    index_max_conv : str
        this will gives you the highest convertion index.

    '''
    sotune1=np.array(data['time'])
    sotune2=np.array(data['temp'])
    sotune3=np.array(data['pressure'])
    sotune4=np.array(data['conv'])
    if app.upper()=='PLOT_TEMP':
       plt.plot(sotune2, sotune1,color='black')
       plt.title('temprature over time')
       plt.xlabel('time(second)')
       plt.ylabel('temprature(celsious) ')
       plt.grid()
       plt.show()
     
    elif app.upper()=='PLOT_PRESSURE':
       plt.plot(sotune3, sotune1,color='red')
       plt.title('pressure over time')
       plt.xlabel('time(second)')
       plt.ylabel('pressure(bar) ')
       plt.grid()
       plt.show()
       
    elif app.upper()=='PLOT_CONVERSION':
          plt.plot(sotune4, sotune1,color='blue')
          plt.title('conversion over time')
          plt.xlabel('time(second)')
          plt.ylabel('conversion(%) ')
          plt.grid()
          plt.show()

    elif app.upper()=='MAXIMUM CONVERSION':
        maxstress=sotune4.max()
        maxstrain=sotune2.max()
        print('maximum of tempreture is ' , maxstrain )
        print('maximum of conversion is ' , maxstress )
        index_max_conv=np.argmax(sotune4) 
        print('The tempreture in maximum conversion is ' , sotune2[index_max_conv] ,
              'and the preesure is ' ,sotune3[index_max_conv]   )
        return index_max_conv
    
    else:
        raise TypeError ('The datas or application is not entered correctly.')
    return sotune1
    return sotune2
    return sotune3  
    return sotune4        

Conversion(data_dasti,'maximum conversion')
Conversion(data_dasti,'plot_pressure')



#import section
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#====================================================================================
#====================================================================================

#########------>>>>>>> Step0 Cleanin Data <<<<<<<<<<<<<----------------

#$$ we can make a function to import a csv file and make sure the data pass the step0 (or cleaning data) 


def Import_Data(File_Directory=None):
    
    
    if File_Directory is None:
        raise TypeError('Please enter the file directory to open!')
        
        
        
    try:
        data=pd.read_csv(File_Directory)
        
        
#for examle -->data=pd.read_csv('C:\\Users\\Parsis.Co\\Desktop\\CV.csv')
         
         
#1--->File not find  (by 'FileNotFoundError')       
    except FileNotFoundError:
        raise FileNotFoundError('Unfortunately, the desired file <',File_Directory,'>is not available ')
        
        
        
#2---> File is  empty (by 'pd.errors.EmptyDataError')
    except pd.errors.EmptyDataError:
         raise ValueError('The file: ',data, ' is empty. please a correct file.')
         
   
#3--->Format is wrong  (by 'pd.errors.ParserError')     
    except pd.errors.ParserError:
             raise ValueError('The file format is not valid, please import a <csv> format file')
             
             
 #4--->remove empty cells     
    if data.isnull().values.any():
        print('Empty cell founded and removed ')
        data.dropna(inplace=True)
        
        
 #5--->turn object to numeric form for both columns    
    for x in data['P']:
        if data['P'].dtype=='object':
            print('object element founded in potential column and converted to numeric form')
            data['P']=pd.to_numeric(data['P'])
            
    for y in data['C']:
        if data['C'].dtype=='object':
            print('object element founded in current density column and converted to numeric')
            data['P']=pd.to_numeric(data['P'])
            
            
 #6--->remove duplicated data         
    if data.duplicated().any():
        print('Duolicated elemets in rows founded and removed ')
        data=data.drop_duplicates()
        
        
    return data
        
            
            
            
 #====================================================================================
 #====================================================================================

 ##font definitions for title and labels

title_style={'family':'times new roman',
             'color':'red',
             'size':28}
label_style={'family':'times new roman',
              'color':'black',
             'size':18}

 #====================================================================================
 #====================================================================================           
            
 
#########------>>>>>>> Main function for cyclic voltammetry investigations <<<<<<<<<<<<<----------------            
            

def CV(data=None,Application=None, is_reversible=None):
    '''
    

    Parameters
    ----------
    data : DataFrame
        Enter your .csv format file as pd.DataFrame.
        Data is the potential vs. current density expoted from the potentiostate device to study electrochemical of your fabricated electrode.
        To use this function correctly, the potential column should named as 'P' and your current density column should named as 'C'
        
    Application : str
        Please enter the application you want to do for your data, including: plot, maxcurrent, peakpotential, charge, and Diffusion-coefficient.
    is_reversible : bool, optional
        If your reaction id reversible please enter 'True', and if the reaction is irreversible enter 'False' to calculate the Diffusion coefficient for your reaction.

    Application=plot (type=plot)--> Plot a cyclic voltammogram (current density vs. potential)
    Application= maxcurrent (type=float)--> the function returns the peak of current that the current density id maximum.
    Application= peakpotential (type=float)--> the function returns the potential which attributed to maxmum current or peak current.
    Application=charge (type=float)--> The function returns the charge corresponding to the integration of the plot.
    Application=Diffusion_coefficient (type=float)--> The function returns the value according the andles-Sevcik aquation depends on reversiblly of irreversibly of reaction.
    

    '''
    
#--> ensure the defined-file name (data) is in calling section
    
    if data is None:
        raise ValueError('Please enter the file name and make sure you import your file as DataFrame')
        
#--> ensure the user enter the application in calling section
        
    if Application is None:
        raise ValueError('Please enter an application for your CV data')
        
    
    if Application=='plot':
        potential=np.array(data['P'])
        current_density=np.array(data['C'])
        plt.plot(potential,current_density,c='c', linewidth=3)
        plt.title('Cyclic Voltammetry',fontdict=title_style)
        plt.xlabel('Potential (V) vs. Ag/AgCl',fontdict=label_style)
        plt.ylabel('Current Density (mA/cm^2)',fontdict=label_style)
        plt.grid(color='gray',linewidth=0.5)
        plt.show()
        
        
    elif Application=='maxcurrent':
        maxcurrent=data['C'].max()
        print('The maximum current density value is: ', maxcurrent, 'mA/cm^2')
        return maxcurrent
        
    
    elif Application=='peakpotential':
        maxcurrent=data['C'].max()
        maxcurrentindex=data['C'].idxmax()
        peakpotetial=data.loc[maxcurrentindex,'P']
        print('The peak potential corresponding to the maximum current is: ',peakpotetial)
        return peakpotetial
        
    elif Application=='charge':
        potential=np.array(data['P'])
        current_density=np.array(data['C'])
        charge=np.trapezoid(current_density,potential)
        print('The charge for your CV is= ', charge, 'c')
        return charge
        
    
    
        
    elif Application=='Diffusion_Coefficient':  #(Peak_Current,A,C,Scan_Rate,n,is_reversible):
        maxcurrent=data['C'].max()
        Area=0.0123
        C=0.2
        Scan_Rate=50
        n=1
        
        if is_reversible is None:
            raise ValueError ('Please enter <True> if the reaction is reversible, else enter <False>')
            
        if is_reversible:
            Diffusion_Coefficient=(maxcurrent/((2.65*10**5)*(n**(3/2))*Area*C*(Scan_Rate**(1/2))))**2
            
        else:
            Diffusion_Coefficient=(maxcurrent/(0.446*(n**(3/2))*Area*C*(Scan_Rate**(1/2))))**2
            
            
        print('The value of Diffusion coefficient for your electrochemical reaction is:',Diffusion_Coefficient )
        return Diffusion_Coefficient

#import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#---------------------- Example 1 ------------------------------------

def Product_Of_Vectors(data,app):
    
    '''
    do bordar (two vectors) be onvane voroodi migirad, va anha ra baham zarb mikonad.
    
    Parameters
    ----------
    data : float
        yek matrix ba do sotoon ast ke sotoone aval, bordar avval, va digari bordar dovvom.
    app : str
          plot agar bashe, rasm mikoneh, ertebat do bordar ro ba ham va agar
          calculate bood zarb dakheli bordar avval ba tranahadeh bordar dovvom
          ra mohasebeh mikonad.

    Returns
    -------
    output : float
        Production of vectors.

    '''
    x=np.array(data['vector1'])
    y=np.array(data['vector2'])
    
    if app=='plot':
        plt.plot(x,y,marker='<',ms=10,mfc='tab:orange',mec='c',c='tab:blue',linewidth=2)
        plt.title('nemoodar')
        plt.xlabel('vector1')
        plt.ylabel('vector2')
        plt.show()
    elif app=='calculate':
        y1=np.transpose(y)
        output=np.dot(x,y1)
        return output


a=np.random.uniform(-5,10,size=100).reshape(-1,1)
b=np.random.uniform(0,10,size=100).reshape(-1,1)
c=np.concatenate([a,b],axis=1)
data=pd.DataFrame(c,columns=['vector1','vector2'])
Product_Of_Vectors(data,'plot') 
Product_Of_Vectors(data,'calculate')   

#-------------------- Example 2 ---------------------------------
def Echelon_Matrix(data,ap):
    '''
    یک ماتریس 90در 200 داریم که مثلا 200 مولفه را در مورد سلامت نود نفر جمع آوری کردیم
    و حالا میخواهیم ببینیم که کدوم یکی از این مولفه ها اضافی هستند و اطلاعات خوبی 
    از سلامتی افراد نمیده. هدفمون اینه که ماتریس را  به صورت سطری پلکانی
    در آوریم یعنی بالا مثلثیش کنیم و مولفه های مهم رو پیدا کنیم. 
    Parameters
    ----------
    data : float
        This is a matrix with more than 50 rows and columns.
    ap : str
          plot agar bashe, rasm mikoneh, ertebat sotoone 10 , 50 ro ba ham va agar
          up bood matrix ra be soorat satri pelekani mikonad.

    Returns
    -------
    A : float
        This is a triangular matrix.

    '''
    x=np.array(data[:,10])
    y=np.array(data[:,50])
    if ap=='plot':
        plt.plot(x,y,marker='4',ms=20,mec='y',c='tab:purple',ls='--')
        plt.title('nemoodare ertebat ghand khoon va vazn')
        plt.xlabel('vazn')
        plt.ylabel('ghand khoon')
        plt.show()
    elif ap=='up':
        A=np.triu(data)
        return A
   
aa=np.random.randint(0,100,(90,200)) 
Echelon_Matrix(aa,'plot')
Echelon_Matrix(aa,'up')





def MyFunction(data,application):
     if application=='plot':
         x=data['DT']
         #print(x)
         y=data['Gold']
         plt.plot(x,y)
         plt.xlabel('Date')
         plt.xticks(rotation=90)
         plt.ylabel('Gold Price per OUNS')
         plt.title('Gold Price base on USD')
         plt.show()  
     elif application=='operation':
         Converting_Gold_Price_to_Toman(data)
         
         
'''
Converting_Gold_Price_to_Toman

input:
data pandas dataframe 
 

output:
data pandas dataframe
'''
def Converting_Gold_Price_to_Toman(ddf):
     a=np.array(ddf['Gold'])
     ddf['Doller'] = pd.to_numeric(ddf['Doller'], errors='coerce')
     b=np.array(ddf['Doller'])
     c=np.multiply(a,b)
     c=np.divide(c,31.1035)
     c=np.multiply(c,1000)
     ddf['Gold price in Toman']=c
     datadata.rename(columns={'Date':'DT', 'GLD': 'Gold', 'DLR': 'Doller'}, inplace=True)
     print('\n')
     print(ddf[['DT','Gold price in Toman']])
     return(ddf)
    


#load data into a DataFrame object from csv file:
datadata=pd.read_csv('Gold_price.csv')
datadata.rename(columns={'Date':'DT', 'GLD': 'Gold', 'DLR': 'Doller'}, inplace=True)
command=0
while(int(command)!=1 and int(command)!=2):
     command=input('\nWelcom to converting Gold Price in first three monthes of 2017: \nwhich action would you like to do?\n1-for drawing a plot press 1\n2-for doing action press 2  ')
if(int(command)==1):
   MyFunction(datadata,'plot')
elif(int(command)==2):
   MyFunction(datadata,'operation')



def Fatigue_Test_Analysis(data,application):
    '''
    

    Parameters
    ----------
    data : data is the exel file with the two columns (stress_amplitude column and number_of_cycles column)
    application : plot , max stress amplitude , fatigue strength , fatigue life , stress in one cycle , Sa , 
                    fatigue limit , std stress , std cycles
    

    Returns
    -------
    plot: S-N plot
    fatigue strength: استحکام خستگی 
    fatigue life: عمر خستگی
    stress in one cycle: Basquin's equation to define B constant.The B is the value of the stress at one cycle.
    Sa: max stress amplitude in cycle.
    fatigue limit: The time of cycle that stress no change.
    std stress: انحراف معیار تنش‌ها
    std cycles: انحراف معیارا سیکل‌ها

    '''
    stress_amplitude = np.array(data["stress_amplitude"])
    number_of_cycles = np.array(data["number_of_cycles"])
    
 
    if application=="plot":
        title_font={"color":"black",
                 'family':'Merriweather',
                 'size': 20     }
        
        xy_label_font={"color":"Magenta",
                 'family':'Merriweather',
                 'size': 12     }
        
        plt.plot(number_of_cycles,stress_amplitude,marker='o',c='c',mec='k',mfc='r',label='S-N Curve',linewidth=3)
        plt.title('S-N Curve (Fatigue Test)',fontdict=title_font,pad=13)
        plt.xscale('log')  
        plt.xlabel('Number of Cycles to Failure (N)',fontdict=xy_label_font)
        plt.ylabel('Stress Amplitude (MPa)',fontdict=xy_label_font)
        plt.grid()
        plt.legend()
        plt.show()
    
    if application=='max stress amplitude' :
        max_stress_amplitude=np.max(stress_amplitude)
        return max_stress_amplitude
    
    if application=='fatigue strength' :
        fatigue_strengt=np.mean(stress_amplitude)   
        return fatigue_strengt
    
    if application=="fatigue life" :
        fatigue_life = np.mean(number_of_cycles)
        return fatigue_life
   
    if application=='stress in one cycle' :
        n=np.log10(number_of_cycles[0])-np.log10(number_of_cycles[1])
        m=np.log10(stress_amplitude[1])-np.log10(stress_amplitude[0])
        slope=n/m
        stress_in_one_cycle=number_of_cycles[0]*2*((stress_amplitude[0])**slope)
        return stress_in_one_cycle
    
    if application=='Sa' :
        i=np.where(number_of_cycles==1)
        Sa=stress_amplitude[i]
        return Sa
   
    if application=='fatigue limit' :
        repetative_index = None
        for i in range(0,len(stress_amplitude)-1) :
            if stress_amplitude[i]==stress_amplitude[i+1] and repetative_index==None :
                repetative_index=i
            elif stress_amplitude[i]!=stress_amplitude[i+1] :
                repetative_index=None
        fatigue_limit=number_of_cycles[repetative_index]
        return fatigue_limit                
        
    if application=='std stress' :
        std_stress = np.std(stress_amplitude)
        return std_stress
     
    if application=='std cycles' :
        std_cycles = np.std(number_of_cycles)
        return std_cycles
    



import openpyxl
import pandas as pd
import numpy as np
import math 
from matplotlib import pyplot as plt

Data = pd.read_excel(r"C:\Users\Bartar\OneDrive\Desktop\project\BTC_Price(2022)2.xlsx" )

Data.info()

Data.dropna(inplace=0)
#a = pd.DataFrame(Data)
#a['Price'] = a['Price'].astype(float)
#a['Data'] = a['Data'].astype(str)

#for i in a.index:
    # data.loc[5,['Temp']]=None
 #   if a.loc[i, ['Price']] < 0:
  #      a = a.drop(i)
   # if DAta.

def Price_BTC(Data,application):
    '''
    
    Parameters
    ----------
    Data : int
        This function needs an excel with two columns (Price and Date)
    application : str
        These are some ways you can use this function  :
            plot
            Max_Price
            Max_Date
            Min_Price
            Min_Price
            

    Returns
    -------
    Plot and int
        plot of price per date
        Date of Maximum and minimum Price

    '''
    Price = np.array(Data['Price'])
    Date = np.array(Data['Date'])
    APP = application.upper()
    
    if APP == 'PLOT' :
        plt.plot(Date,Price,ms=3,mfc='b',c='k',linewidth=3)
        plt.title('Price of BTC')
        plt.ylabel('(US)Price')
        plt.xlabel('Date(AD)')
        plt.show()
        
    elif APP =='MAX_PRICE' :
        return Date[np.argmax(Price)], Price.max()
    
    elif APP == 'MAX_DATE' :
        return Date[-1], Price[-1]
    
    elif APP == 'MIN_PRICE' :
        return Date[np.argmin(Price)] , Price.min()
    
    elif APP == 'MIN_DATE':
        return Date[0],Price[0]
    
print (Price_BTC(Data,'plot'))












