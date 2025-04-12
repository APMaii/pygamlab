# PyGamMLab

**PyGAMLab** is a scientific Python library developed for researchers, engineers, and students who need access to fundamental constants, conversion tools, engineering formulas, and data analysis utilities. The package is designed with simplicity, clarity, and usability in mind.

---

## 📌 Overview

**PyGAMLab** stands for *Python GAMLAb tools*, a collection of scientific tools and functions developed at the **GAMLab** under the supervision of **Prof. Nader** at **Amirkabir University of Technology (AUT)**.

- **Author:** Ali Pilehvar Meibody  
- **Supervisor:** Prof. Nader  
- **Affiliation:** GAMLab, Amirkabir University of Technology (AUT)

---

## 📦 Modules

PyGAMLab is composed of **four core modules**, each focused on a specific area of scientific computation:

### 🔹 `constants.py`
This module includes a comprehensive set of scientific constants used in physics, chemistry, and engineering.

Examples:
- Planck's constant
- Boltzmann constant
- Speed of light
- Universal gas constant
- And many more...

---

### 🔹 `convertors.py`
Contains unit conversion functions that follow the format:  
`FirstUnit_To_SecondUnit()`

Examples:
- `Kelvin_To_Celsius(k)`
- `Pascal_To_Bar(p)`
- `Meter_To_Foot(m)`
- ...and many more standard conversions used in science and engineering.

---

### 🔹 `functions.py`
This module provides a wide collection of **scientific formulas and functional tools** commonly used in engineering applications.

Examples:
- Thermodynamics equations
- Mechanical stress and strain calculations
- Fluid dynamics formulas
- General utility functions

---

### 🔹 `data_analysis.py`
Provides tools for working with data, either from a **file path** or directly from a **DataFrame**.

Features include:
- Reading and preprocessing datasets
- Performing scientific calculations
- Creating visualizations (e.g., line plots, scatter plots, histograms)

---

## 🚀 Installation

To install PyGAMLab via pip (after uploading it to PyPI):

```bash
pip install pygamlab
```


---

## 📖 Usage Example

```python
from pygamlab import Kelvin_To_Celsius, universal_gas_constant
from pygamlab import calculate_stress
from pygamlab import plot_from_csv

# Convert temperature
print(Kelvin_To_Celsius(300))

# Use a constant
print(universal_gas_constant)

# Use an engineering function
stress = calculate_stress(force=500, area=0.01)
print(stress)

# Analyze a dataset
plot_from_csv("data/experiment.csv", column_x="Time", column_y="Temperature")
```

