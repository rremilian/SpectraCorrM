import os
import re
import numpy as np
import scipy as sp
import csv
from rich.table import Table
from rich.console import Console
from matplotlib import pyplot as plt

"""
SpectraCorrM- A Python module for the analysis of experimental and theoretical IR and Raman spectra 
@github: https://github.com/rremilian
"""

class LogFile:

    def __init__(self, location, raman=False):
        self.location = location
        self.raman = raman
        if os.name == "nt":
            self.name = self.location.split("\\")[-1]
        elif os.name == "posix":
            self.name = self.location.split("/")[-1]
        self.freqlist = []
        self.intlist = []
    
    def display(self):
        table = Table(title=self.name)
        table.add_column("Mode", justify="right")
        table.add_column("Frequency (cm^-1)", justify="center")
        table.add_column("Intensity", justify="left")
        
        for index, values in enumerate(zip(self.freqlist, self.intlist)):
            table.add_row(str(index+1), str(values[0]), str(values[1]))
            
        console = Console()
        console.print(table)
    
    def spectrum(self,emin, emax, step, sigma, scale=1.00, mode="lorentz"):
        return TheoreticalSpectrum(self.freqlist, self.intlist, emin, emax, step, sigma, scale, mode)

class LogFileGau(LogFile):
        
    def __init__(self, location, raman):
        LogFile.__init__(self, location, raman)
        
        # Extract frequencies, IR intensities and Raman activities from the log file
        with open(self.location, "r") as file:
            lines = file.read().split("\n")
            for line in lines:
                if "Frequencies ---" in line:
                    line = line.split()
                    [self.freqlist.append(float(freq)) for freq in line[2:7]]
                if self.raman == True and "Raman Activities ---" in line:
                    line = line.split()
                    [self.intlist.append(float(ramanint)) for ramanint in line[3:8]]
                elif self.raman == False and "IR Intensities ---" in line:
                    line = line.split()
                    [self.intlist.append(float(irint)) for irint in line[3:8]]
            file.close()

        # Raise exception if len(freqlist) != len(ramanintens)
        if raman == True and (len(self.freqlist) != len(self.intlist)):
            raise Exception("The number of frequencies is not equal to the number of raman activities")
        if raman == False and (len(self.freqlist) != len(self.intlist)):
            raise Exception("The number of frequencies is not equal to the number of IR intensities")

class LogFileORCA(LogFile):
    
    def __init__(self, location, raman):
        LogFile.__init__(self, location, raman)
        
        # Extract frequencies, IR intensities and Raman activities from the log file
        with open(self.location, "r") as file:
            lines = file.read().split("\n")
            write_raman = False
            write_ir = False
            for line in lines:
                if raman == True and "RAMAN SPECTRUM" in line:
                    write_raman = True
                elif raman == False and "IR SPECTRUM" in line:
                    write_ir = True
                if write_raman and re.search("\d:",line):
                    line = line.split()
                    self.freqlist.append(float(line[1]))
                    self.intlist.append(float(line[2]))
                elif "The first" in line:
                    write_raman = False
                if write_ir and re.search("\d:", line):
                    line = line.split()
                    self.freqlist.append(float(line[1]))
                    self.intlist.append(float(line[3]))
                elif "* The epsilon" in line:
                    write_ir = False
            file.close()
        
        # Raise exception if len(freqlist) != len(ramanintens)
        if raman == True and (len(self.freqlist) != len(self.intlist)):
            raise Exception("The number of frequencies is not equal to the number of raman activities")
        if raman == False and (len(self.freqlist) != len(self.intlist)):
            raise Exception("The number of frequencies is not equal to the number of IR intensities")

class Spectrum:
    
    def __init__(self, frequencies, intensities):
        self.frequencies = frequencies
        self.intensities = intensities
    
    @classmethod
    def from_csv(cls, location):
        frequencies = []
        intensities = []
        with open(location, "r") as file:
            csv_file = csv.reader(file, delimiter=",")
            for line in csv_file:
                if re.search("\d",line[0]) and re.search("\d",line[1]):
                    frequencies.append(float(line[0]))
                    intensities.append(float(line[1]))
                else:
                    pass
        file.close()
        
        if frequencies[2] < frequencies[1]:
            frequencies = frequencies[::-1]
            intensities = intensities[::-1]

        return cls(frequencies, intensities)
    
    def normalize(self):
        max_intensity = self.intensities.max()
        self.intensities /= max_intensity

    def plot(self):
        plt.plot(self.frequencies, self.intensities)
        plt.show()
    
    def export_csv(self, location):
        header = ["Frequency", "Intensity"]
        with open(location, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(header)
            intensities_round = [round(intensity, 4) for intensity in self.intensities]
            writer.writerows(zip(self.frequencies, intensities_round))
    
    def export_png(self, location, xlabel="Frequency", ylabel="Intensity", dpi=300):
        plt.plot(self.frequencies, self.intensities)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(location, bbox_inches='tight', dpi=dpi)
    
    def interpolate(self, thspectrum):
        interpolated_intensities = np.interp(thspectrum.frequencies, self.frequencies, self.intensities)
        return interpolated_intensities
    
    def compare(self, thspectrum, mode="pearson"):
        interpolated_intensities = Spectrum.interpolate(self, thspectrum)
        if mode == "pearson":
            r = sp.stats.pearsonr(interpolated_intensities, thspectrum.intensities)[0]
        elif mode == "spearman":
            r = sp.stats.spearmanr(interpolated_intensities, thspectrum.intensities)[0]
        return r
            
class TheoreticalSpectrum(Spectrum):

    def __init__(self, freqlist, intlist, emin, emax, step, sigma, scale, raman, mode="lorentz"):
        self.freqlist = freqlist
        self.intlist = intlist
        
        # Check if the variables have the correct type
        for param in [emin, emax, step, sigma, scale]:
            if isinstance(param, str):
                raise TypeError("emax, emin, step, sigma and scale must have type int or float")
            
        self.emin = emin
        self.emax = emax

        # Check if emin > emax
        if self.emin > self.emax:
            raise Exception("The maximum frequency is lower than the minimum frequency")
        
        self.sigma = sigma
        self.step = step
        self.nstep = int(round((self.emax - self.emin)/self.step)+1)
        self.raman = raman

        self.scale = scale
        if self.scale != 1.00:
            self.freqlist_scaled = [i*self.scale for i in self.freqlist]
        else:
            self.freqlist_scaled = self.freqlist

        # Writing frequencies values to a list
        self.frequencies = []
        for i in range(self.nstep):
            self.frequencies.append(self.emin + i * self.step)

        # Define distributions
        self.mode = mode
        if self.mode == "gauss":
            gauss = lambda params, x: (1.0/(params[1]*np.sqrt(2*np.pi))*np.exp(-0.5*((x-params[0])/params[1])**2))
        elif self.mode == "lorentz":
            lorentz  = lambda params, x: (1.0/np.pi*params[1])/((x-params[0])**2+params[1]**2)
        else:
            raise Exception("The mode is incorrect. It should be gauss or lorentz")

        # Generate the spectrum
        temp = np.empty((self.nstep, len(self.freqlist)))
        for i in range(len(self.freqlist_scaled)):
            params = [self.freqlist_scaled[i],self.sigma]
            if self.mode == "lorentz":
                temp[:,i] = lorentz(params, np.asarray(self.frequencies, dtype=float)) * self.intlist[i]
            elif self.mode == "gauss":
                temp[:,i] = gauss(params, np.asarray(self.frequencies, dtype=float)) * self.intlist[i]
            self.intensities = np.sum(temp, axis=1, dtype=float)
