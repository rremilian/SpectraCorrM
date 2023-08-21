from setuptools import setup, find_packages

setup(
    name = "Spectra Analyze",
    version = "1.0.0",
    url = "https://github.com/rremilian/spectra_analyze",
    author = "rremilian",
    description = "A Python module for the analysis of experimental and theoretical IR and Raman spectra ",
    packages=find_packages(),
    install_requires= ["numpy >= 1.21.5", "matplotlib >= 3.5.1", "scipy >= 1.7.3", "rich >= 13.3.5"],
)
