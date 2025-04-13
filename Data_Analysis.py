'''
Data_Analysis.py :









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
import numpy as np
import pandas as pd
#import matplotlib as plt
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns 
from scipy import stats
from scipy.integrate import solve_bvp



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
         average_pressure=statistics.mean(pressure)
         return average_pressure
        
    elif which=='average volume':
         average_volume=statistics.mean(volume)
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
        draw : Drawing a Plot
        calculate : Calculate the Consupmtion Energie in [mWs]
        please say which work we do ( 1 or 2).

    '''
        
    mydf=pd.read_excel(input_file)
    
    
    if which=='draw':
        #get the data on each columns
        
        A=mydf['Angle[°]']
        
        Energie =mydf['Energie']
       
        #plotting data
        plt.plot(A, Energie,color = 'green')
        plt.title('Energie of OTT Motor 185000 Cycle')
        plt.xlabel('Angle[°]')
        plt.ylabel('Consupmtion Energie')
        plt.show()
    
    
    if which=='calculate':
        mydf=pd.DataFrame(mydf,columns=['Angle[°]','Power[mW]','Time for a Cycle','Energie'])
        
        summ = mydf['Energie'].sum()                          # The amount of Energie for a half Cycle of Duty in mWs
       
        summ =( summ * 2)/1000                                # The amount of Consumption Energie for a Dutycycle in Ws
        
        return summ
        







def Stress_Strain1(df,operation,L0=90,D0=9):
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
       It tells the function to whether PLOT the curve or find the UTS valu. 
       
     L0: initial length of the sample
     D0: initial diameter of the sample

    Returns
    -------
    The Stress-Strain curve or the amount of UTS
    
    '''
    
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







def Stress_Strain2(input_file,which,count):
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
    
    
def Stress_Strain3(input_data, action):
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

def Stress_Strain4(file_path, D0, L0):
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




def Stress_Strain5(input_data, action):
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



def aerospace_analysis (CSV,which):
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
    ages_x = np.arange(18, 49)

    # get income data
    python_income = python_data['Income']
    js_income = js_data['Income']

    # get max income data
    max_income_python = python_data['Income'].max()
    max_income_js = js_data['Income'].max()



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
    



def DMTA_Test(Data2,Operator,Sample_Name):
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




def Wear_Rate(w,work,S,F):
    '''
    w : dataframe contains 2 columns 
    bayad shamele 2 soton bashad ke dar yek soton vazn nemoone ghabl az test andazegiri shavad ,
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
    









def Polarization(b,work):
    '''
    d: dadehaye azmayeshgah shmele ghegaliye jaryan va potansiel hastand.
    '''
    
    font_1={'family': 'serif',
            'color': 'b',
            'size':14}

    font_2={'family': 'serif',
            'color': 'k',
            'size':12}
    
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
            max intensity
            scatter plot
            line graph
            line graph fill between
    Returns
    -------
    max intensity and degree of max intensity 
    or
    plot

    '''
    
    font_title = { 'family' : 'serif','color' :'black' , 'size' : 15 }
    font_x = {'family' :  'serif','color' :   'black', 'size' : 10 }
    font_y = { 'family' : 'serif','color' :   'black', 'size' : 10}    


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
    font_title = { 'family' : 'serif','color' :'black' , 'size' : 15 }
    font_x = {'family' :  'serif','color' :   'black', 'size' : 10 }
    font_y = { 'family' : 'serif','color' :   'black', 'size' : 10}   
    
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
         img.show()




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


   




def Brain_Data_example():
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




'''

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

'''

#________________________________Flux Calculation______________________________
'''

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


'''







def SI_Calculation(f_loc,P,PC,Density=1):

    

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







def Porosity(e_loc,Density=1):

    

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




def Tortuosity(e_loc):

    

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







def Pore_Size(g_loc,A,P,Vis=8.9*1e-4):

    

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








def Membrane_Rejection(Experimental_Datas):
    Dox_Size=1.68*10**-9 #m
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









def XRD_Analysis2(file,which,peak=0):
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






def aerospace2 (CSV,which):
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
                   
                   
   
                    
    
    
    
   
def Energie2(input_file,which):
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
    


 

def fatigue(f_loc,which):
    new_data=pd.read_csv(f_loc)
    
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




def XRD4(data,application):
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
        



def Stress_Strain6(data,application):
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
    title_style={'family':'times new roman',
                 'color':'red',
                 'size':28}
    label_style={'family':'times new roman',
                  'color':'black',
                 'size':18}
    
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
    







