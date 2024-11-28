import numpy as np

class ProfileFunction:
    """
    Initialize through
    param::list [alpha,beta,gamma] values described in our papers
    param::float between 0,1 (not the 0 itself)
    """
    def __init__(self, parameters, x):
        # Calculate the function at initialization
        self.funcAtPoint = parameters[0] * np.power(1 - x, 3) * np.log(1 / x) \
                         + parameters[1] * np.power(1 - x, 3) \
                         + parameters[2] * x * np.power(1 - x, 2)

    def __call__(self):
        # Recompute the function with new parameters and x
        return self.funcAtPoint
    

class deltaProfileFunction:
    """
    Initiliaze through the following params and call () to get the values.
    param::string analysis Set e.g. "Set11"
    param::string GPD type e.g. "Ht"
    param::string flavor e.g. "uv"
    param::float between 0,1 (not the 0 itself)
    """
    def __init__(self, AnalysisSet, gpdType, flavor,x):
        self.AnalysisSet = AnalysisSet
        self.gpdType = gpdType
        self.flavor = flavor 
        self.x = x

        if "Set11" == self.AnalysisSet:
            if "H" == self.gpdType:
                if "uv" == self.flavor:
                    self.uncertainty =self.__get_delta_Set11_H_uv__(self.x)
                elif "dv" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_H_dv__(self.x)
                elif "ubar" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_H_ubar__(self.x)
                elif "dbar" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_H_dbar__(self.x)
            elif "Ht" == self.gpdType:
                if "uv" == self.flavor:
                    self.uncertainty =self.__get_delta_Set11_Ht_uv__(self.x)
                elif "dv" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_Ht_dv__(self.x)
                elif "ubar" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_Ht_ubar__(self.x)
                elif "dbar" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_Ht_dbar__(self.x)
            elif "E" == self.gpdType:
                if "uv" == self.flavor:
                    self.uncertainty =self.__get_delta_Set11_E_uv__(self.x)
                elif "dv" == self.flavor:
                    self.uncertainty =  self.__get_delta_Set11_E_dv__(self.x)



    def __call__(self):
        # Recompute the function with new parameters and x
        return self.uncertainty
    

    ##################


############################### GPD H SET 11 ###############################

    def __get_delta_Set11_H_uv__(self,x):
        term1 = 0.000853492921079 * (1. - x)**6
        term2 = -0.001895848197352 * (1. - x)**5 * x
        term3 = 0.001262260111969 * (1. - x)**4 * x**2
        term4 = -0.000421501341966 * (1. - x)**6 * np.log(1. / x)
        term5 = 0.000439598050862 * (1. - x)**5 * x * np.log(1. / x)
        term6 = 0.000055937850344 * (1. - x)**6 * (np.log(1. / x)**2)
        return np.sqrt(0. + term1 + term2 + term3 + term4 + term5 + term6)
    


    def __get_delta_Set11_H_dv__(self, x):
        term1 = 0.019268853290179 * (1 - x)**6
        term2 = -0.112728975221076 * (1 - x)**5 * x
        term3 = 0.202003814039234 * (1 - x)**4 * x**2
        term4 = -0.006391431559782 * (1 - x)**6 * np.log(1 / x)
        term5 = 0.016575965450532 * (1 - x)**5 * x * np.log(1 / x)
        term6 = 0.000580804456069 * (1 - x)**6 * (np.log(1 / x))**2

        result = np.sqrt(0 + term1 + term2 + term3 + term4 + term5 + term6)
        return result

    def __get_delta_Set11_H_ubar__(self,x):
        term1 = 0.000853492921079 * (1 - x)**6
        term2 = -0.001895848197352 * (1 - x)**5 * x
        term3 = 0.001262260111969 * (1 - x)**4 * x**2

        result = np.sqrt(0 + term1 + term2 + term3)
        return result

    def __get_delta_Set11_H_dbar__(self, x):
        term1 = 0.019268853290179 * (1 - x)**6
        term2 = -0.112728975221076 * (1 - x)**5 * x
        term3 = 0.202003814039234 * (1 - x)**4 * x**2

        result = np.sqrt(0 + term1 + term2 + term3)
        return result

############################### GPD HT SET 11 ###############################


    def __get_delta_Set11_Ht_uv__(self, x):
        term1 = 0.046690116436714 * (1 - x)**6
        term2 = -0.400339308478034 * (1 - x)**5 * x
        term3 = 0.973836212063479 * (1 - x)**4 * x**2
        term4 = -0.000584662433724 * (1 - x)**6 * np.log(1 / x)
        term5 = 0.00076340954987 * (1 - x)**5 * x * np.log(1 / x)
        term6 = 0.000055937850344 * (1 - x)**6 * (np.log(1 / x))**2

        result = np.sqrt(term1 + term2 + term3 + term4 + term5 + term6)
        return result

    def __get_delta_Set11_Ht_dv__(self, x):
        term1 = 0.007572274268172 * (1 - x)**6
        term2 = -0.04323940902673 * (1 - x)**5 * x
        term3 = 0.135807387348546 * (1 - x)**4 * x**2
        term4 = -0.004007433991774 * (1 - x)**6 * np.log(1 / x)
        term5 = 0.007781531430436 * (1 - x)**5 * x * np.log(1 / x)
        term6 = 0.000580804456069 * (1 - x)**6 * (np.log(1 / x))**2

        result = np.sqrt(term1 + term2 + term3 + term4 + term5 + term6)
        return result

    def __get_delta_Set11_Ht_ubar__(self, x):
        term1 = 0.046690116436714 * (1 - x)**6
        term2 = -0.400339308478034 * (1 - x)**5 * x
        term3 = 0.973836212063479 * (1 - x)**4 * x**2
        term4 = -0.000584662433724 * (1 - x)**6 * np.log(1 / x)
        term5 = 0.00076340954987 * (1 - x)**5 * x * np.log(1 / x)
        term6 = 0.000055937850344 * (1 - x)**6 * (np.log(1 / x))**2

        result = np.sqrt(term1 + term2 + term3 + term4 + term5 + term6)
        return result

    def __get_delta_Set11_Ht_dbar__(self, x):
        term1 = 0.007572274268172 * (1 - x)**6
        term2 = -0.04323940902673 * (1 - x)**5 * x
        term3 = 0.135807387348546 * (1 - x)**4 * x**2
        term4 = -0.004007433991774 * (1 - x)**6 * np.log(1 / x)
        term5 = 0.007781531430436 * (1 - x)**5 * x * np.log(1 / x)
        term6 = 0.000580804456069 * (1 - x)**6 * (np.log(1 / x))**2

        result = np.sqrt(term1 + term2 + term3 + term4 + term5 + term6)
        return result






############################### GPD E SET 11 ###############################

    def __get_delta_Set11_E_uv__(self, x):
        term1 = 0.01333874899347 * (1. - x)**6
        term2 = -0.012154794598044 * (1. - x)**5 * x
        term3 = 0.029208166459806 * (1. - x)**4 * x**2
        term4 = -0.015938384805344 * (1. - x)**6 * np.log(1. / x)
        term5 = -0.010947625571654 * (1. - x)**5 * x * np.log(1. / x)
        term6 = 0.014809707405028 * (1. - x)**6 * (np.log(1. / x)**2)
        return np.sqrt(0. + term1 + term2 + term3 + term4 + term5 + term6)


    def __get_delta_Set11_E_dv__(self, x):
        term1 = 0.013824883857035 * (1 - x)**6
        term2 = -0.113441561962804 * (1 - x)**5 * x
        term3 = 0.448898192980077 * (1 - x)**4 * x**2
        term4 = -0.006385664770452 * (1 - x)**6 * np.log(1 / x)
        term5 = 0.030653416150456 * (1 - x)**5 * x * np.log(1 / x)
        term6 = 0.001217329821894 * (1 - x)**6 * (np.log(1 / x))**2

        result = np.sqrt(0 + term1 + term2 + term3 + term4 + term5 + term6)
        return result