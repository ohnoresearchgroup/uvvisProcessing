import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapezoid as trapz


class SolFlSpectrum():
    def __init__(self,fullpath):
        self.fullpath = fullpath

        
        data = pd.read_csv(fullpath,skiprows=12,sep='\t',names =['wl','counts','extra'])
        data = data.drop(columns='extra')
        data = data.dropna()


        self.data = data
        self.bg = np.mean(self.data['counts'][0:10])
        self.data['flcorr']=self.data['counts']-self.bg
        self.data['flnorm']=self.data['flcorr']/np.max(self.data['flcorr'])
        
        plt.figure()
        plt.plot(self.data['wl'],self.data['flcorr'])
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Fluorescence [arb.u.]')
        
        print(data)
        
        
    def saveSpectrum(self,outputpath):
        self.data.to_csv(outputpath + "PROCESSED_" + self.name,index = False)


    def integrate_spectrum(self,wlrange):
        wl1 = wlrange[0]
        wl2 = wlrange[1]

        idx1 = np.argmin(np.abs(self.data['wl'] - wl1))
        idx2 = np.argmin(np.abs(self.data['wl'] - wl2))
        
        return trapz(self.data['counts'][idx1:idx2], self.data['wl'][idx1:idx2])
        