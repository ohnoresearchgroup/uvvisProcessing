import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class AbsSpectrum():
    def __init__(self,fullpath,horiba=False):
        self.fullpath = fullpath

        if horiba is False:
            data = pd.read_csv(fullpath,skiprows=2, names =['wl','abs','extra'] )
            data = data.drop(columns='extra')
            data = data.dropna()
        else:
            data = pd.read_csv(fullpath,skiprows=12,sep='\t',names =['wl','abs','extra'])
            data = data.drop(columns='extra')
            data = data.dropna()


        self.data = data
        self.bg = np.mean(self.data['abs'][0:10])
        self.data['abscorr']=self.data['abs']-self.bg
        self.data['absnorm']=self.data['abscorr']/np.max(self.data['abscorr'])
        
        plt.figure()
        plt.plot(self.data['wl'],self.data['abscorr'])
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Absorbance [a.u.]')
        
        print(data)
        
        
    def saveSpectrum(self,outputpath):
        self.data.to_csv(outputpath + "PROCESSED_" + self.name,index = False)