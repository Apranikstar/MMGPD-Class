import numpy as np
from scipy.integrate import quad
from .profileFuncClass import ProfileFunction 
from .csvParserClass import getProfileFunctionParameters
class Observables:

    def __init__(self, mmgpdDOTgpdAnalysis):
        self.__gpdAnalysis = mmgpdDOTgpdAnalysis
        self.m_p =  0.93827208943
        self.m_n =  0.93956542052
        self.m_p2 = np.power(self.m_p,2)
        self.m_n2 = np.power(self.m_n,2)


    def d1_with_sea(self, t,xi):
        return np.divide(self.__term_d1_with_sea__(t,xi),xi**2) * np.divide(5,4)

    def d1_without_sea(self, t,xi):
        return np.divide(self.__term_d1_without_sea__(t,xi),xi**2) * np.divide(5,4)

#########################
    def r2mass_p_w_A0(self, D0):
        
        result = (6 * self.__diffA0__() - (3 * np.divide(D0, 2 * self.m_p2))) * np.divide(1, self.__A0__())
        return result * (0.1973 ** 2)
    
    def r2mass_p_wo_A0(self, D0):
        
        result = (6 * self.__diffA0__() - (3 * np.divide(D0, 2 * self.m_p2))) 
        return result * (0.1973 ** 2)


##############################################################################################################################
  ########################################################Subroutines########################################################
##############################################################################################################################

                ########################### d1 Subroutines ###########################
    def __integrand_terms_d1_w_sea__(self, x,t,xi):
        return (-self.__gpdAnalysis.xGPD("Set11", "H" , "uv" , x , t)   -2 * self.__gpdAnalysis.xGPD("Set11", "H" , "ubar" , x , t)
                - self.__gpdAnalysis.xGPD("Set11", "H" , "dv" , x , t)  - 2 * self.__gpdAnalysis.xGPD("Set11", "H" , "dbar" , x , t) 
                + self.__gpdAnalysis.xGPDxi("Set11", "H" , "uv" , x , t,xi)   +2 * self.__gpdAnalysis.xGPDxi("Set11", "H" , "ubar" , x , t,xi)
                + self.__gpdAnalysis.xGPDxi("Set11", "H" , "dv" , x , t,xi)  + 2 * self.__gpdAnalysis.xGPDxi("Set11", "H" , "dbar" , x , t,xi) 
               )

    def __term_d1_with_sea__(self,t,xi):
        return  quad (self.__integrand_terms_d1_w_sea__ , 0 , 1 , args = (t,xi) , limit = 250)[0]
    

    def __integrand_terms_d1_wout_sea__(self, x,t,xi):
        return (-self.__gpdAnalysis.xGPD("Set11", "H" , "uv" , x , t)   
                - self.__gpdAnalysis.xGPD("Set11", "H" , "dv" , x , t)  
                + self.__gpdAnalysis.xGPDxi("Set11", "H" , "uv" , x , t,xi)   
                + self.__gpdAnalysis.xGPDxi("Set11", "H" , "dv" , x , t,xi)  
               )
    
    def __term_d1_without_sea__(self,t,xi):
        return  quad (self.__integrand_terms_d1_wout_sea__ , 0 , 1 , args = (t,xi) , limit = 250)[0]
    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################

                ########################### r2 mass Subroutines ###########################

    
    def __A0__(self):
        def __A0integrand__(x,t):
            return (self.__gpdAnalysis.xGPD("Set11", "H" , "uv" , x , t)   +2 * self.__gpdAnalysis.xGPD("Set11", "H" , "ubar" , x , t)
                    + self.__gpdAnalysis.xGPD("Set11", "H" , "dv" , x , t)  + 2 * self.__gpdAnalysis.xGPD("Set11", "H" , "dbar" , x , t) )
        t = 0
        return quad (__A0integrand__ , 0 , 1 , args = (t,) , limit = 250)[0]
    

    def __diffA0__(self):
        def __integranddiffA__(x):
            profFuncParameters={}
            flavor = ["uv","dv","ubar", "dbar"]
            for items in flavor:
                profFuncParameters[items] = getProfileFunctionParameters("HGAG23", "H", "Set11")(items)
            return (      ProfileFunction(profFuncParameters["uv"],x)()    *self.__gpdAnalysis.xGPD("Set11", "H" , "uv" , x , 0)
                +     ProfileFunction(profFuncParameters["dv"],x)()    *self.__gpdAnalysis.xGPD("Set11", "H" , "dv" , x , 0) 
                + 2 * ProfileFunction(profFuncParameters["ubar"],x)()  *self.__gpdAnalysis.xGPD("Set11", "H" , "ubar" , x , 0)
                + 2 * ProfileFunction(profFuncParameters["dbar"],x)()  *self.__gpdAnalysis.xGPD("Set11", "H" , "dbar" , x , 0)
               )
        return quad (__integranddiffA__ , 0 , 1, limit = 250)[0]

    ##############################################################################################################################
    ##############################################################################################################################
    ##############################################################################################################################
