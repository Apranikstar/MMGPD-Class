import numpy as np
from scipy.integrate import quad
from .xPDFClass import xPDF


class EMObservables:
    
    def __init__(self, pdfSET):
        self.pdf = xPDF(pdfSET)
        self.__chargeUV = np.divide(2,3)
        self.__chargeDV = np.divide(-1,3)
        self.__chargeSV = np.divide(-1,3)
        self.__m_p  = 0.93827
        self.__m_n = 0.93956

########################################################
    def GM(self, t, Q2, #x,
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv,
           ID,  
           ):
        
        F1uv = self.__chargeUV * quad(self.__H__, 1e-9, 1, args=(Q2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv"),limit = 200)[0] 
        F1dv = self.__chargeDV * quad(self.__H__, 1e-9, 1, args=(Q2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv"),limit = 200)[0] 
        F1sv = self.__chargeSV * quad(self.__H__, 1e-9, 1, args=(Q2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv"),limit = 200)[0] 

        F2uv = quad(self.__E__, 1e-9, 1, args=(t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv"),limit = 200)[0]
        F2dv = quad(self.__E__, 1e-9, 1, args=(t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv"),limit = 200)[0]
        F2sv = quad(self.__E__, 1e-9, 1, args=(t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv,"sv"),limit = 200)[0]

        if 1 == ID:
            F1p = self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv 
            F2p = self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv 
            return  F1p + F2p

        if 2 == ID:
            F1n = self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv 
            F2n = self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F1sv 
            return F1n + F2n


########################################################
    def GE(self, t, Q2, #x,
           H_aprime_uv, H_aprime_dv, H_aprime_sv,
           H_B_uv, H_B_dv, H_B_sv,
           H_A_uv, H_A_dv, H_A_sv,
           E_aprime_uv, E_aprime_dv, E_aprime_sv,
           E_B_uv, E_B_dv, E_B_sv,
           E_A_uv, E_A_dv, E_A_sv,
           E_alpha_uv, E_alpha_dv, E_alpha_sv,
           E_beta_uv, E_beta_dv, E_beta_sv,
           E_gamma_uv, E_gamma_dv,E_gamma_sv  ,
           ID
           ):
        

        F1uv = self.__chargeUV * quad(self.__H__, 1e-9, 1, args=(Q2,t, H_aprime_uv, H_B_uv, H_A_uv,"uv"),limit = 200)[0] 
        F1dv = self.__chargeDV * quad(self.__H__, 1e-9, 1, args=(Q2,t,H_aprime_dv, H_B_dv, H_A_dv, "dv"),limit = 200)[0] 
        F1sv = self.__chargeSV * quad(self.__H__, 1e-9, 1, args=(Q2,t, H_aprime_sv, H_B_sv, H_A_sv,"sv"),limit = 200)[0] 

        F2uv = quad(self.__E__, 1e-9, 1, args=(t, E_aprime_uv, E_B_uv, E_A_uv, E_alpha_uv, E_beta_uv, E_gamma_uv, "uv"),limit = 200)[0]
        F2dv = quad(self.__E__, 1e-9, 1, args=(t, E_aprime_dv, E_B_dv, E_A_dv, E_alpha_dv, E_beta_dv, E_gamma_dv, "dv"),limit = 200)[0]
        F2sv = quad(self.__E__, 1e-9, 1, args=(t, E_aprime_sv, E_B_sv, E_A_sv, E_alpha_sv, E_beta_sv, E_gamma_sv,"sv"),limit = 200)[0]

        if 1 == ID:
            F1p = self.__chargeUV * F1uv + self.__chargeDV * F1dv + self.__chargeSV * F1sv 
            F2p = self.__chargeUV * F2uv + self.__chargeDV * F2dv + self.__chargeSV * F2sv 
            return F1p + np.divide(t,4 * self.__m_p**2) * F2p

        if 2 == ID:
            F1n = self.__chargeDV * F1uv + self.__chargeUV * F1dv + self.__chargeSV * F1sv 
            F2n = self.__chargeDV * F2uv + self.__chargeUV * F2dv + self.__chargeSV * F1sv 
            return F1n + np.divide(t,4* self.__m_n **2) * F2n



##################################  Subroutines ##################################

    def __f__(self, x, aprime, B, A):
        return aprime *  np.power(1- x,3) * np.log( np.divide(1,x) ) + B * np.power(1-x,3) + A * x * np.power(1-x,2)
    
    def __H__(self, x,Q2,t,aprime,B,A,flavor):
        return self.pdf.xPDFCentVal(flavor,x,Q2)/x * np.exp(t * self.__f__(x,aprime,B,A))
    
    def __E__(self, x, t, aprime, B, A, alpha, beta, gamma, flavor):
        return self.__pdfEHandler__(flavor, x,alpha,beta,gamma) * np.exp(t * self.__f__(x,aprime,B,A))
    
    def __pdfEHandler__(self, flavor, x,alpha,beta,gamma):
        def integrand(x, alpha, beta, gamma):
            return (np.power(x, -alpha) * np.power(1-x, beta) * (1 + gamma * np.sqrt(x)))
        k = {
        "uv": 1.67,
        "dv": -2.03,
        "sv": 0
            
        }
        N = np.divide(1, quad(integrand, 1e-9, 1, args=(alpha, beta, gamma))[0])
        return k.get(flavor) * N * np.power(x,-alpha) * np.power(1-x,beta) * (1+ gamma * np.sqrt(x))
