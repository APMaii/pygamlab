from setuptools import setup, find_packages




setup(
    name="pygamlab",
    version="2.0.0",
    author="Ali Pilehvar Meibody",
    author_email="alipilehvar1999@gmail.com",
    description="A scientific Python library with constants, unit convertors, formulas, and data analysis tools.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/APMaii/pygamlab/tree/main",
    license="MIT",  # Correct field (ensure no 'license-file' field)
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "numpy", "pandas", "scipy", "matplotlib", "seaborn", "scikit-learn",
        "ase", "plotly", "PyQt5", "jarvis", "mp_api", "aflow"
    ],
    include_package_data=True,
    package_data={
        "pygamlab.ai_core": ["gam_models/*.gam_ai"],
    },
)
