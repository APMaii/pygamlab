'''
Constants.py ==> 

This module contains scientific constants intended for use throughout the package.
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

'''
import math

# ==============================================================================
# Physical and Scientific Constants
# ==============================================================================

# A
A = 0.413  # lattice constant in nanometers
alpha = 5.67e-8  # Stefan-Boltzmann constant (W/m²/K⁴)
Avogadro_Number = 6.022e23  # Avogadro's number (mol⁻¹)

# B
Boltzman_constany = 1.381e-23  # Boltzmann constant (J/K)
Boltzmann_Constant = 1.380649e-23  # More precise value (J/K)

# C
c = 2.998e8  # Speed of light in vacuum (m/s)
C = 299792458  # Speed of light (duplicate)
Conductivity_P3HT = 2.4  # Electrical conductivity of P3HT
Conductivity_PLN = 10  # Electrical conductivity of PLN
Conductivity_PPY = 105  # Electrical conductivity of polypyrrole
Conductivity_pT = 33.7  # Electrical conductivity of polythiophene

# D
d = 1.128  # Statistical constant for n = 2
D = 2.3e-5  # Oxygen diffusion coefficient in water (µm²/s)
density_of_Al = 2.7  # g/cm³
density_of_Cu = 8.96  # g/cm³
density_of_Fe = 7.87  # g/cm³

# E
Earth_Accel = 9.8  # Acceleration due to gravity on Earth (m/s²)
Electron_Charge = 1.6e-19  # Elementary charge (C)
e = 2.718281828459045  # Euler's number
Eutectic_Percent = 4.3  # Carbon % in eutectic reaction
Eutectic_T = 1148  # Eutectic temperature (°C)
Eutectoid_persent = 0.76  # Carbon % in eutectoid reaction
Eutectoid_T = 727  # Eutectoid temperature (°C)

# F
F = 96485  # Faraday constant (C/mol)
Faraday_Constant = 96485  # Duplicate (C/mol)
Fe_Density = 7.87  # Density of iron
Fe_Tm_Alpha = 910  # Melting point of alpha-phase iron (°C)
Fe_Tm_Delta = 1539  # Melting point of delta-phase iron (°C)
Fe_Tm_Gama = 1495  # Melting point of gamma-phase iron (°C)

# G
G = 6.674e-11  # Gravitational constant (m³/kg/s²)
G_Mol_Ba = 137.33  # Molar mass of barium
G_Mol_O = 16.00  # Molar mass of oxygen
G_Mol_Si = 28.09  # Molar mass of silicon
G_Mol_Ti = 47.87  # Molar mass of titanium
G_Mol_Zr = 91.22  # Molar mass of zirconium
g = 9.81  # Standard gravity (m/s²)
g0 = 9.81  # Duplicate of g

# H
h = 6.62607015e-34  # Planck’s constant (J·s)
H = 6.634e-34  # Less precise duplicate

# I
# (None explicitly declared)

# K
K = 1.38e-23  # Boltzmann constant (J/K)
k = 1.381e-23  # Duplicate Boltzmann constant
k = 534.6  # Also reused, possibly context-specific
k = 1.96  # Z value for standard normal distribution
Kw = 1e-14  # Ion product of water at 25°C

# L
Latent_heat = 1.16e9  # Latent heat (unit not specified)

# M
max_C_inSteel = 2.11  # Max carbon % in steel before it becomes cast iron
melting_point_of_Al = 660  # °C
melting_point_of_Cu = 1085  # °C
melting_point_of_Fe = 1538  # °C

# N
N_A = 6.022e23  # Avogadro’s constant
Na = 6.022e23  # Duplicate
N_a = 6.022e23  # Duplicate
N_t = 3e33  # Possibly total atoms or particles (context needed)

# P
P_0 = 101325  # Standard atmospheric pressure (Pa)
phi = 1.618  # Golden ratio
Pi = 3.141592653589793  # Pi
π = 3.14  # Duplicate with symbol

# Q
q = 0.7  # Possibly heat transfer or general coefficient

# R
R = 8.314  # Ideal gas constant (J/mol·K)
R_Cal = 1.987  # Gas constant in cal/(mol·K)
R_LA = 0.08205  # Ideal gas constant in (L·atm)/(mol·K)

# S
S = 28.34  # Standard entropy for solid aluminum (J/mol·K)
Speed_Of_Light = 2.99e8  # Another duplicate of c

# T
T_m = 1064  # Melting temperature in °C
ta = 6.28  # Approximate tau
Tau = 2 * Pi  # Tau (2π)
thermal_conductivity_coefficient_of_Al = 237  # W/m·K
thermal_conductivity_coefficient_of_Cu = 385  # W/m·K
thermal_conductivity_coefficient_of_Fe = 221  # W/m·K
Tm = 14.025  # Melting point of hydrogen in K
t_stu = 0.9277  # t-Student value

# V
vacuum_permeability = 1  # Simplified value

# X, Z
Xd = 0.95
Xw = 0.05
Zf = 0.25
zeta = 1.202  # Riemann zeta function at 3 (used in physics)
