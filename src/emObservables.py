import numpy as np
from scipy.integrate import quad
from .xPDFClass import xPDF
from scipy.special import gamma as gamma_f   # gamma function



class EMObservables:
    
    def __init__(self, pdfSET,mu2, mp, mn , mup , mun):
        self.__pdf = xPDF(pdfSET)
        self.__chargeUV = np.divide(2,3)
        self.__chargeDV = np.divide(-1,3)
        self.__chargeSV = np.divide(-1,3)
        self.__mup = mup
        self.__mun = mun
        self.__m_p  = mp
        self.__m_n = mn
        self.__mu2 = mu2
        self.__alpha_qed = 1.0 / 137.035999084
        #def NORM(x,flavor,mu2):
        #    return self.__pdf.xPDFCentVal(flavor,x,mu2)/x
        #Nuv = quad(NORM, 1e-9 , 1, args=("uv",self.__mu2),limit = 300)[0]
        #Ndv = quad(NORM, 1e-9 , 1, args=("dv",self.__mu2), limit = 300)[0]
        ##Nsv = quad(NORM, 1e-9 , 1, args=("sv",self.__mu2))[0]
        Nuv = quad(lambda x: self.__pdf.xPDFCentVal("uv", x, self.__mu2) / x, 1e-9, 1, limit=150)[0]
        Ndv = quad(lambda x: self.__pdf.xPDFCentVal("dv", x, self.__mu2) / x, 1e-9, 1, limit=150)[0]
        
        self.__NormDict = {
            "uv" : np.divide(2 , Nuv),
            "dv" : np.divide(1 , Ndv),
            "sv" : 1
        }
        
########################################################
    def G_D(self, t):
        return  np.array([(  1 - t[i] / 0.71 )**(-2) for i in range(len(t))])
        
        
    def GM(self, ID, t, #mu2, 
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv,
           Ks  
           ):
        self.__ks = Ks
        ### Flavor Form Factors
        def F1(x,ID,t):
            F1uv =  self.__H__(x, self.__mu2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv")
            F1dv =  self.__H__(x, self.__mu2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv")
            F1sv =  self.__H__(x, self.__mu2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv
            if 2 == ID:
                return self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv
                
        def F2(x,ID,t):
            F2uv =  self.__E__(x, t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv")
            F2dv =  self.__E__(x, t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv")
            F2sv = self.__E__(x, t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv
            if 2 == ID:
                return self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F2sv

        def F1F2GM(x,ID,t):
            return F1(x,ID,t) + F2(x,ID,t)
        if 1 == ID:
            return  np.array([quad(F1F2GM,1e-9,1, args=(1,t[i]), limit = 250)[0] for i in range(len(t))])

        if 2 == ID:
            return np.array( [quad(F1F2GM,1e-9,1, args=(2,t[i]) ,limit = 250 )[0] for i in range(len(t))])
############################################################
    def GE(self, ID, t, #mu2, 
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv,
           Ks  
           ):
        self.__ks = Ks
        ### Flavor Form Factors
        def F1(x,ID,t):
            F1uv =  self.__H__(x, self.__mu2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv")
            F1dv =  self.__H__(x, self.__mu2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv")
            F1sv =  self.__H__(x, self.__mu2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv
            if 2 == ID:
                return self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv
            
                
        def F2(x,ID,t):
            F2uv =  self.__E__(x, t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv")
            F2dv =  self.__E__(x, t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv")
            F2sv = self.__E__(x, t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv
            if 2 == ID:
                return self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F2sv
                
        def F1F2GE(x, ID, m,t):
            return F1(x,ID,t) +  np.divide(t,4 * m**2) * F2(x,ID,t)

        if 1 == ID:
            return np.array([quad(F1F2GE,1e-9,1, args=( 1,self.__m_p,t[i]), limit = 250)[0] for i in range(len(t))])

        if 2 == ID:
            return  np.array([quad(F1F2GE,1e-9,1,  args=( 2,self.__m_n,t[i]) ,limit = 250 )[0] for i in range(len(t))] )  
############################################################
    def RatioGEGM(self, ID, t, #mu2, 
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv,
           Ks  
           ):
        self.__ks = Ks
        ### Flavor Form Factors
        def F1(x,ID,t):
            F1uv =  self.__H__(x, self.__mu2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv")
            F1dv =  self.__H__(x, self.__mu2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv")
            F1sv =  self.__H__(x, self.__mu2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv
            if 2 == ID:
                return self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv
            
                
        def F2(x,ID,t):
            F2uv =  self.__E__(x, t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv")
            F2dv =  self.__E__(x, t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv")
            F2sv = self.__E__(x, t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv
            if 2 == ID:
                return self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F2sv
                
        def F1F2GE(x, ID, m,t):
            return (F1(x,ID,t) +  np.divide(t,4 * m**2) * F2(x,ID,t))
        def F1F2GM(x,ID,t):
            return F1(x,ID,t) + F2(x,ID,t)
        
        if 1 == ID:
            result =  (np.array([quad(F1F2GE,1e-9,1, args=( 1,self.__m_p,t[i]), limit = 250)[0] for i in range(len(t))]))/np.array([quad(F1F2GM,1e-9,1, args=( 1,t[i]), limit = 250)[0] for i in range(len(t))  ])
            return  np.array(result)

        if 2 == ID:
            result =  (np.array([quad(F1F2GE,1e-9,1, args=( 2,self.__m_n,t[i]), limit = 250)[0]for i in range(len(t))]))/np.array([quad(F1F2GM,1e-9,1, args=( 2,t[i]), limit = 250)[0] for i in range(len(t))  ])
            return  np.array(result)

##################################  PRad SigR ##################################
    def SigR(self, ID, t, episilon, #mu2, 
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv,
           Ks,
           ):
        self.__ks = Ks
        ### Flavor Form Factors
        def F1(x,ID,t):
            F1uv =  self.__H__(x, self.__mu2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv")
            F1dv =  self.__H__(x, self.__mu2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv")
            F1sv =  self.__H__(x, self.__mu2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv
            if 2 == ID:
                return self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv
            
                
        def F2(x,ID,t):
            F2uv =  self.__E__(x, t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv")
            F2dv =  self.__E__(x, t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv")
            F2sv =  self.__E__(x, t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv, "sv")
            if 1 == ID:
                return self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv
            if 2 == ID:
                return self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F2sv
                
        def F1F2GE(x, ID, m,t):
            return (F1(x,ID,t) +  np.divide(t,4 * m**2) * F2(x,ID,t))
        def F1F2GM(x,ID,t):
            return F1(x,ID,t) + F2(x,ID,t)
        
        if 1 == ID:
            #tau = np.divide(-t , 4 * np.power(self.__m_p,2)) 
            return   (np.array([episilon[i] * np.power(quad(F1F2GE,1e-9,1, args=( 1,self.__m_p,t[i]), limit = 250)[0],2) for i in range(len(t))]))+np.array([ np.divide(-t[i] , 4 * np.power(self.__m_p,2)) * np.power(quad(F1F2GM,1e-9,1, args=( 1,t[i]), limit = 250)[0],2) for i in range(len(t))  ])
        if 2 == ID:
            #tau = np.divide(-t , 4 * np.power(self.__m_n,2)) 
            return   (np.array([episilon[i] * np.power(quad(F1F2GE,1e-9,1, args=( 2,self.__m_n,t[i]), limit = 250)[0],2) for i in range(len(t))]))+np.array([ np.divide(-t[i] , 4 * np.power(self.__m_n,2)) * np.power(quad(F1F2GM,1e-9,1, args=( 2,t[i]), limit = 250)[0],2) for i in range(len(t))  ])
##################################  PRad SigTot ##################################
    def SigTot(self, ID, t, #mu2, 
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv,
           Ks,
           theta,El  
           ):
        self.__ks = Ks
        ### Flavor Form Factors
        def F1(x,ID,t):
            F1uv =  self.__H__(x, self.__mu2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv")
            F1dv =  self.__H__(x, self.__mu2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv")
            F1sv =  self.__H__(x, self.__mu2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv
            if 2 == ID:
                return self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv
            
                
        def F2(x,ID,t):
            F2uv =  self.__E__(x, t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv")
            F2dv =  self.__E__(x, t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv")
            F2sv = self.__E__(x, t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv,"sv")
            if 1 == ID:
                return self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv
            if 2 == ID:
                return self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F2sv
                
        def F1F2GE(x, ID, m,t):
            return (F1(x,ID,t) +  np.divide(t,4 * m**2) * F2(x,ID,t))
        def F1F2GM(x,ID,t):
            return F1(x,ID,t) + F2(x,ID,t)


        if 1 == ID:
            tau = np.divide(self.__mu2 , 4 * np.power(self.__m_p,2))
            theta_rad = theta * np.pi / 180.0 
            episilon = 1.0 / (1 + 2 * (1 + tau) * np.tan(theta_rad / 2.0)**2)
            denominator = np.divide(1, episilon * (1 + tau)) 
            El_gev = El / 1000.0 
            epl_gev = self.__mu2 / (4 * El_gev * np.sin(theta_rad / 2.0)**2)
            conv = 3.893793e-1 
            sig_mott = (
            self.__alpha_qed**2 / (4 * El_gev**2 * np.sin(theta_rad / 2.0)**4)
            * epl_gev / El_gev * np.cos(theta_rad / 2.0)**2
            )
            dMott = conv * sig_mott
            return  dMott * denominator * (np.array([episilon * quad(F1F2GE,1e-9,1, args=( 1,self.__m_p,t[i]), limit = 250)[0] for i in range(len(t))]))+np.array([tau * quad(F1F2GM,1e-9,1, args=( 1,t[i]), limit = 250)[0] for i in range(len(t))  ])

        if 2 == ID:
            tau = np.divide(self.__mu2 , 4 * np.power(self.__m_n,2))
            theta_rad = theta * np.pi / 180.0 
            episilon = 1.0 / (1 + 2 * (1 + tau) * np.tan(theta_rad / 2.0)**2)
            denominator = np.divide(1, episilon * (1 + tau))
            El_gev = El / 1000.0 
            epl_gev = self.__mu2 / (4 * El_gev * np.sin(theta_rad / 2.0)**2)
            conv = 3.893793e-1 
            sig_mott = (
            self.__alpha_qed**2 / (4 * El_gev**2 * np.sin(theta_rad / 2.0)**4)
            * epl_gev / El_gev * np.cos(theta_rad / 2.0)**2
            ) 
            dMott = sig_mott * conv
            return  dMott * denominator * (np.array([quad(F1F2GE,1e-9,1, args=( 2,self.__m_n,t[i]), limit = 250)[0]for i in range(len(t))]))+np.array([quad(F1F2GM,1e-9,1, args=( 2,t[i]), limit = 250)[0] for i in range(len(t))  ])
              

##################################  Subroutines ##################################

    def __f__(self, x, aprime, B, A):
        return aprime *  np.power(1- x,3) * np.log( np.divide(1,x) ) + B * np.power(1-x,3) + A * x * np.power(1-x,2)
    
    def __H__(self, x,mu2,t,aprime,B,A,flavor):
        return  self.__NormDict.get(flavor) * self.__pdf.xPDFCentVal(flavor,x,mu2)/x * np.exp(t * self.__f__(x,aprime,B,A)) 
    
    def __E__(self, x, t, aprime, B, A, alpha, beta, gamma, flavor):
        return self.__pdfEHandler__(flavor, x,alpha,beta,gamma) * np.exp(t * self.__f__(x,aprime,B,A))
    
    def __pdfEHandler__(self, flavor, x,alpha,beta,gamma):

        # kp= mup-1  ,  kn= mun , 
        #   if (Vks==0.d0)ku= 2*kp + kn , kd= kp + 2*kn
        # else:  ku= 2*kp + kn + Vks  , kd= kp + 2*kn + Vks
         
        if 0 == self.__ks:

            k  = {
            "uv": 2 * (self.__mup - 1) + self.__mun,
            "dv": (self.__mup - 1) + 2 * self.__mun,
            "sv": 0
            }
        elif 0 != self.__ks:
            k = {
            "uv":  2 * (self.__mup - 1) + self.__mun + self.__ks,
            "dv": (self.__mup - 1) + 2 * self.__mun + self.__ks,
            "sv": self.__ks
            }
        N = np.divide(1,(gamma_f(1.+beta) * gamma_f(1.0-alpha)/gamma_f(2.-alpha+beta))  +
             (gamma*gamma_f(1.+beta)*gamma_f(1.5-alpha)/gamma_f(2.5-alpha+beta))   )
        
        return k.get(flavor) * N * np.power(x,-alpha) * np.power(1-x,beta) * (1+ gamma * np.sqrt(x))
