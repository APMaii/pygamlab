from setuptools import setup, find_packages




setup(
    name="PyGamLab",
    version="0.2.1",
    author="Ali Pilehvar Meibody",
    author_email="your_email@example.com",
    description="A scientific Python library with constants, unit convertors, formulas, and data analysis tools.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/APMaii/pygamlab/tree/main",
    license="MIT",  # Correct field (ensure no 'license-file' field)
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "numpy", "pandas", "scipy", "matplotlib", "seaborn", "scikit-learn"
    ],
    #include_package_data=True,
)
