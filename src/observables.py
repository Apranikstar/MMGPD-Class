import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

class Observables:

    def __init__(self, mmgpdDOTgpdAnalysis):
        self.gpdAnalysis = mmgpdDOTgpdAnalysis
        

    def d1(self, t,xi):
        return np.divide(self.__term_d1__(t,xi),xi**2) * np.divide(5,4)

            
    def __integrand_terms_d1__(self, x,t,xi):
        return (-self.gpdAnalysis.xGPD("Set11", "H" , "uv" , x , t)   -2 * self.gpdAnalysis.xGPD("Set11", "H" , "ubar" , x , t)
           - self.gpdAnalysis.xGPD("Set11", "H" , "dv" , x , t)  - 2 * self.gpdAnalysis.xGPD("Set11", "H" , "dbar" , x , t) 
           + self.gpdAnalysis.xGPDxi("Set11", "H" , "uv" , x , t,xi)   +2 * self.gpdAnalysis.xGPDxi("Set11", "H" , "ubar" , x , t,xi)
           + self.gpdAnalysis.xGPDxi("Set11", "H" , "dv" , x , t,xi)  + 2 * self.gpdAnalysis.xGPDxi("Set11", "H" , "dbar" , x , t,xi) )

    def __term_d1__(self,t,xi):
        return  quad (self.__integrand_terms_d1__ , 0 , 1 , args = (t,xi) , limit = 160)[0]
