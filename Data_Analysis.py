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
import cmath
import random
import numpy as np
import pandas as pd
import matplotlib as plt





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









