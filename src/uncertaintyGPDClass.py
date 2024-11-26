import numpy as np

class UncertaintyGPD:
    
    def __init__(self, analysisSet, gpdType, flavor, x, t):
        self.analysisSet = analysisSet 
        self.gpdType =  gpdType
        self.flavor = flavor
        self.x = x
        self.t = t

        if "Set11" == self.analysisSet:
            if "H" == self.gpdType:
                if "uv" == self.flavor:
                    self.uncertainty = self.__SigHuvSet11(self.x, self.t)
                elif "dv" == self.flavor:
                    self.uncertainty = self.__SigHdvSet11(self.x, self.t)
                elif "ubar" == self.flavor:
                    self.uncertainty = self.__SigHubarSet11(self.x, self.t)
                elif "dbar" == self.flavor:
                    self.uncertainty = self.__SigHdbarSet11(self.x, self.t)
            elif "Ht" == self.gpdType:
                if "uv" == self.flavor:
                    self.uncertainty = self.__SigHtuvSet11(self.x, self.t)
                elif "dv" == self.flavor:
                    self.uncertainty = self.__SigHtdvSet11(self.x, self.t)
                elif "ubar" == self.flavor:
                    self.uncertainty = self.__SigHtubarSet11(self.x, self.t)
                elif "dbar" == self.flavor:
                    self.uncertainty = self.__SigHtdbarSet11(self.x, self.t)
            elif "E" == self.gpdType:
                if "uv" == self.flavor:
                    self.uncertainty = self.__SigEuvSet11(self.x, self.t)
                elif "dv" == self.flavor:
                    self.uncertainty = self.__SigEdvSet11(self.x, self.t)


#################### SET 11 Subroutines #####################


############# GPD H #############
    def __SigHuvSet11(self, x,t):
        coeff1 = 1.07041282750141
        coeff2 = 0.806536120991285
        coeff3 = 0.659296317994651

        # Common terms
        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        # Full expression
        Huv = np.sqrt(
        0.0
        + 0.000853492921079 * t**2 * one_minus_x**6 * exp_term
        - 0.001895848197352 * t**2 * one_minus_x**5 * x * exp_term
        + 0.001262260111969 * t**2 * one_minus_x**4 * x**2 * exp_term
        - 0.000421501341966 * t**2 * one_minus_x**6 * exp_term * log_term
        + 0.000439598050862 * t**2 * one_minus_x**5 * x * exp_term * log_term
        + 0.000055937850344 * t**2 * one_minus_x**6 * exp_term * log_term**2)
        return Huv
    
    def __SigHdvSet11(self, x,t):
        coeff1 = 1.29268529185374
        coeff2 = 3.5246852793991
        coeff3 = 0.408689992464628

        # Common terms
        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        # Full expression
        SigHdv = np.sqrt(
        0.0
        + 0.019268853290179 * t**2 * one_minus_x**6 * exp_term
        - 0.112728975221076 * t**2 * one_minus_x**5 * x * exp_term
        + 0.202003814039234 * t**2 * one_minus_x**4 * x**2 * exp_term
        - 0.006391431559782 * t**2 * one_minus_x**6 * exp_term * log_term
        + 0.016575965450532 * t**2 * one_minus_x**5 * x * exp_term * log_term
        + 0.000580804456069 * t**2 * one_minus_x**6 * exp_term * log_term**2)
        return SigHdv
    
    def __SigHubarSet11(self, x,t):
        coeff1 = 1.07041282750141
        coeff2 = 0.806536120991285
        coeff3 = 2.0  # Explicit coefficient for the log term

        # Common terms
        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        # Full expression
        SigHubar = np.sqrt(
        0.0
        + 0.000853492921079 * t**2 * one_minus_x**6 * exp_term
        - 0.001895848197352 * t**2 * one_minus_x**5 * x * exp_term
        + 0.001262260111969 * t**2 * one_minus_x**4 * x**2 * exp_term
        )
        return SigHubar

    def __SigHdbarSet11(self, x,t):
        coeff1 = 1.29268529185374
        coeff2 = 3.5246852793991
        coeff3 = 2.0  # Explicit coefficient for the log term

        # Common terms
        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        # Full expression
        SigHdbar = np.sqrt(
        0.0
        + 0.019268853290179 * t**2 * one_minus_x**6 * exp_term
        - 0.112728975221076 * t**2 * one_minus_x**5 * x * exp_term
        + 0.202003814039234 * t**2 * one_minus_x**4 * x**2 * exp_term
        )
        return SigHdbar


############# GPD HT #############
    def __SigHtuvSet11(self, x,t):
        coeff1 = -0.207814555359072
        coeff2 = 5.4421177732338
        coeff3 = 0.659296317994651

        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        SigHTuv = np.sqrt(
        0.0
        + 0.046690116436714 * t**2 * one_minus_x**6 * exp_term
        - 0.400339308478034 * t**2 * one_minus_x**5 * x * exp_term
        + 0.973836212063479 * t**2 * one_minus_x**4 * x**2 * exp_term
        - 0.000584662433724 * t**2 * one_minus_x**6 * exp_term * log_term
        + 0.00076340954987 * t**2 * one_minus_x**5 * x * exp_term * log_term
        + 0.000055937850344 * t**2 * one_minus_x**6 * exp_term * log_term**2
        )
        return SigHTuv

    def __SigHtdvSet11(self, x,t):
        coeff1 = -1.16651740634566
        coeff2 = 4.89231766427256
        coeff3 = 0.408689992464628

        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        SigHTdv = np.sqrt(
        0.0
        + 0.007572274268172 * t**2 * one_minus_x**6 * exp_term
        - 0.04323940902673 * t**2 * one_minus_x**5 * x * exp_term
        + 0.135807387348546 * t**2 * one_minus_x**4 * x**2 * exp_term
        - 0.004007433991774 * t**2 * one_minus_x**6 * exp_term * log_term
        + 0.007781531430436 * t**2 * one_minus_x**5 * x * exp_term * log_term
        + 0.000580804456069 * t**2 * one_minus_x**6 * exp_term * log_term**2
        )
        return SigHTdv


    def __SigHtubarSet11(self, x,t):
        coeff1 = -0.207814555359072
        coeff2 = 5.4421177732338
        coeff3 = 0.659296317994651

        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        SigHTubar = np.sqrt(
        0.0
        + 0.046690116436714 * t**2 * one_minus_x**6 * exp_term
        - 0.400339308478034 * t**2 * one_minus_x**5 * x * exp_term
        + 0.973836212063479 * t**2 * one_minus_x**4 * x**2 * exp_term
        - 0.000584662433724 * t**2 * one_minus_x**6 * exp_term * log_term
        + 0.00076340954987 * t**2 * one_minus_x**5 * x * exp_term * log_term
        + 0.000055937850344 * t**2 * one_minus_x**6 * exp_term * log_term**2
        )
        return SigHTubar



    def __SigHtdbarSet11(self, x,t):
        coeff1 = -1.16651740634566
        coeff2 = 4.89231766427256
        coeff3 = 0.408689992464628

        one_minus_x = 1.0 - x
        log_term = np.log(1.0 / x)
        exp_term = np.exp(2 * t * (coeff1 * one_minus_x**3 + coeff2 * one_minus_x**2 * x + coeff3 * one_minus_x**3 * log_term))
    
        SigHTdbar = np.sqrt(
        0.0
        + 0.007572274268172 * t**2 * one_minus_x**6 * exp_term
        - 0.04323940902673 * t**2 * one_minus_x**5 * x * exp_term
        + 0.135807387348546 * t**2 * one_minus_x**4 * x**2 * exp_term
        - 0.004007433991774 * t**2 * one_minus_x**6 * exp_term * log_term
        + 0.007781531430436 * t**2 * one_minus_x**5 * x * exp_term * log_term
        + 0.000580804456069 * t**2 * one_minus_x**6 * exp_term * log_term**2
        )
        return SigHTdbar


############# GPD E #############


    def __SigEuvSet11(self, x,t):
        exp_common = np.exp(t * (-0.620521094076626 * (1 - x)**3 +
                             0.966675366165284 * (1 - x)**2 * x +
                             1.025462513758 * (1 - x)**3 * np.log(1 / x)))
        logx = np.log(1 / x)
        log1_minus_x = np.log(1 - x)
    
        term1 = 0.10828783719862625 * t**2 * (1 - x)**14.25299617988472 * (1 - x)**6 * x**1.03615603867022 * \
            np.exp(2 * t * (-0.620521094076626 * (1 - x)**3 + 0.966675366165284 * (1 - x)**2 * x + 
            1.025462513758 * (1 - x)**3 * logx))
    
        term2 = -0.09867615165860649 * t**2 * (1 - x)**14.25299617988472 * (1 - x)**5 * x**2.03615603867022 * \
            np.exp(2 * t * (-0.620521094076626 * (1 - x)**3 + 0.966675366165284 * (1 - x)**2 * x + 
            1.025462513758 * (1 - x)**3 * logx))
    
        term3 = 0.23712037583271445 * t**2 * (1 - x)**14.25299617988472 * (1 - x)**4 * x**3.03615603867022 * \
            np.exp(2 * t * (-0.620521094076626 * (1 - x)**3 + 0.966675366165284 * (1 - x)**2 * x + 
            1.025462513758 * (1 - x)**3 * logx))
    
        term4 = -0.1293924355166352 * t**2 * (1 - x)**14.25299617988472 * (1 - x)**6 * x**1.03615603867022 * \
            np.exp(2 * t * (-0.620521094076626 * (1 - x)**3 + 0.966675366165284 * (1 - x)**2 * x + 
            1.025462513758 * (1 - x)**3 * logx)) * logx
    
        term5 = -0.08887600300411579 * t**2 * (1 - x)**14.25299617988472 * (1 - x)**5 * x**2.03615603867022 * \
            np.exp(2 * t * (-0.620521094076626 * (1 - x)**3 + 0.966675366165284 * (1 - x)**2 * x + 
            1.025462513758 * (1 - x)**3 * logx)) * logx
    
        term6 = 0.1202295046724442 * t**2 * (1 - x)**14.25299617988472 * (1 - x)**6 * x**1.03615603867022 * \
            np.exp(2 * t * (-0.620521094076626 * (1 - x)**3 + 0.966675366165284 * (1 - x)**2 * x + 
            1.025462513758 * (1 - x)**3 * logx)) * logx**2
    
        term7 = -0.227488508657512 * t * (1 - x)**7.12649808994236 * (1 - x)**3 * x**0.51807801933511 * exp_common * \
            (0.18701319953666662 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common + 
             2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * log1_minus_x)
        
        term8 = -0.26815289787482854 * t * (1 - x)**7.12649808994236 * (1 - x)**2 * x**1.51807801933511 * exp_common * \
            (0.18701319953666662 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common + 
             2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * log1_minus_x)
    
        term9 = 0.2148089118499239 * t * (1 - x)**7.12649808994236 * (1 - x)**3 * x**0.51807801933511 * exp_common * logx * \
            (0.18701319953666662 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common + 
             2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * log1_minus_x)
    
        term10 = 0.312379255787054 * (0.18701319953666662 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common + 
                                  2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * 
                                  log1_minus_x)**2
    
        term11 = 0.026598349810500986 * t * (1 - x)**7.12649808994236 * (1 - x)**3 * x**0.51807801933511 * exp_common * \
             (-11.325627623636962 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common - 
              2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * np.log(x))
    
        term12 = 0.043656421341730686 * t * (1 - x)**7.12649808994236 * (1 - x)**2 * x**1.51807801933511 * exp_common * \
             (-11.325627623636962 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common - 
              2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * np.log(x))
    
        term13 = -0.06413311484414962 * t * (1 - x)**7.12649808994236 * (1 - x)**3 * x**0.51807801933511 * exp_common * logx * \
             (-11.325627623636962 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common - 
              2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * np.log(x))
    
        term14 = -0.066580633797924 * (0.18701319953666662 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common + 
                                   2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * 
                                   log1_minus_x) * \
             (-11.325627623636962 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common - 
              2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * np.log(x))
    
        term15 = 0.009137286714088 * (-11.325627623636962 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common - 
                                  2.849261368721917 * (1 - x)**7.12649808994236 * x**0.51807801933511 * exp_common * 
                                  np.log(x))**2

        return np.sqrt(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + 
                   term11 + term12 + term13 + term14 + term15)





    def __SigEdvSet11(self, x,t):
        # Precompute repetitive terms
        term1 = 1 - x
        term2 = np.exp(2 * t * (-1.32488193988708 * term1**3 + 
                            3.78040091664933 * term1**2 * x + 
                            0.743114214316566 * term1**3 * np.log(1. / x)))
        term3 = np.exp(t * (-1.32488193988708 * term1**3 + 
                        3.78040091664933 * term1**2 * x + 
                        0.743114214316566 * term1**3 * np.log(1. / x)))
        log_x = np.log(1. / x)
        log_1_x = np.log(1 - x)
    
        # Components of the equation
        part1 = 0.09668431387536028 * t**2 * term1**35.7773470540264 * term1**6 * x**0.76070102768486 * term2
        part2 = -0.7933534702167954 * t**2 * term1**35.7773470540264 * term1**5 * x**1.76070102768486 * term2
        part3 = 3.139369143132609 * t**2 * term1**35.7773470540264 * term1**4 * x**2.76070102768486 * term2
        part4 = -0.04465814131632227 * t**2 * term1**35.7773470540264 * term1**6 * x**0.76070102768486 * term2 * log_x
        part5 = 0.21437464061838032 * t**2 * term1**35.7773470540264 * term1**5 * x**1.76070102768486 * term2 * log_x
        part6 = 0.008513395107470952 * t**2 * term1**35.7773470540264 * term1**6 * x**0.76070102768486 * term2 * log_x**2
    
        # Nested terms for clarity
        inner_term1 = -0.05412811039177612 * term1**17.8886735270132 * x**0.38035051384243 * term3
        inner_term2 = -2.6445224856167133 * term1**17.8886735270132 * x**0.38035051384243 * term3 * log_1_x
    
        part7 = (0.041501477027317206 * t * term1**17.8886735270132 * term1**3 * x**0.38035051384243 * term3 *
             (inner_term1 + inner_term2))
        part8 = (5.72192832695901 * t * term1**17.8886735270132 * term1**2 * x**1.38035051384243 * term3 *
             (inner_term1 + inner_term2))
        part9 = (0.2264377443796347 * t * term1**17.8886735270132 * term1**3 * x**0.38035051384243 * term3 * log_x *
             (inner_term1 + inner_term2))
        part10 = (8.96776057655549 * (inner_term1 + inner_term2)**2)
        inner_term3 = (14.92358366559217 * term1**17.8886735270132 * x**0.38035051384243 * term3 +
                   2.644522485616713 * term1**17.8886735270132 * x**0.38035051384243 * term3 * np.log(x))
    
        part11 = (-0.0011591799552562803 * t * term1**17.8886735270132 * term1**3 * x**0.38035051384243 * term3 *
              inner_term3)
        part12 = (-0.038741000507643744 * t * term1**17.8886735270132 * term1**2 * x**1.38035051384243 * term3 *
              inner_term3)
        part13 = (0.0007776341779881862 * t * term1**17.8886735270132 * term1**3 * x**0.38035051384243 * term3 * log_x *
              inner_term3)
        part14 = (-0.084434803086666 * (inner_term1 + inner_term2) * inner_term3)
        part15 = (0.000616162279468 * inner_term3**2)
    
        # Combine all parts
        return np.sqrt(part1 + part2 + part3 + part4 + part5 + part6 + 
                   part7 + part8 + part9 + part10 + 
                   part11 + part12 + part13 + part14 + part15)
    


