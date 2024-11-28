import lhapdf
import os
from .csvParserClass import getProfileFunctionParameters
from .profileFuncClass import ProfileFunction
import numpy as np
from scipy.integrate import quad
from uncertainties import ufloat
from .uncertaintyGPDClass import UncertaintyGPD



class GPDAnalysis:
    """
    Class used to access the following:
    xGPD(analysisSet, gpdType, flavor, x, t)
    xGPDwUnc(analysisSet, gpdType, flavor, x, t)
    xGPDxi(analysisSet, gpdType, flavor, x, t,xi)
    """
    """
        Initialize the analysis with a specific Analysis.
        param::string analysisName e.g "HGAG23"
    """
    def __init__(self, analysis_type):
        """
        Initialize the analysis with a specific Analysis.
        param::string analysisName e.g "HGAG23"
        """
        self.name = analysis_type
        self.__analysis_type = analysis_type
        self.Q2 = self.__get_Q2__()
        self.__UPDF = self.__get_analysis_updf__()
        self.__PPDF = self.__get_analysis_ppdf__()
        self.__UGridPDFSet = self.__select_ugrid_pdf_set__()
        self.__PGridPDFSet = self.__select_pgrid_pdf_set__()
        self.__UGridPDFs = self.__get_ugrid_pdfs__()
        self.__PGridPDFs = self.__get_pgrid_pdfs__()
        self.print_analysis_doi()
        self.__flavor_map = {
            "u": 2, "ubar": -2, "uv": (2, -2),
            "d": 1, "dbar": -1, "dv": (1, -1),
            "s": 3, "sbar": -3, "sv": (3, -3),
            "c": 4, "cbar": -4, "cv": (4, -4),
            "b": 5, "bbar": -5, "bv": (5, -5),
            "t": 6, "tbar": -6, "tv": (6, -6),
            "g": 21
            }


    def xGPD(self,analysisSet, gpdType, flavor, x, t):
        """
        param::string::analysisSet e.g. "Set11"
        param::string::gpdType e.g. "Ht"
        param::string::flavor e.g. "dv"
        param::float::x values between 0,1 (not the 0 itself)
        param::float::t Negative values (t=0 returns the forward limit)
        """
        profFuncParameters = getProfileFunctionParameters(self.__analysis_type, gpdType, analysisSet)(flavor)
        profileFunction = ProfileFunction(profFuncParameters,x)() # last () is intentional, don't remove
        if "H" == gpdType:
            pdfFunction = self.__pdfHandler__(flavor,x,self.__UPDF)
        elif "Ht" == gpdType:
            pdfFunction = self.__pdfHandler__(flavor,x,self.__PPDF)
        elif "E" == gpdType:
            pdfFunction = self.__pdfEHandler__(analysisSet, flavor, x)
        return pdfFunction * np.exp(t * profileFunction)
    
    def xGPDwUnc(self,analysisSet, gpdType, flavor, x, t):
        """
        Returns ufloat(xGPD,uncertainty)
        use .n to get the nominal value
        use .s to get the uncertainty value
        param::string::analysisSet e.g. "Set11"
        param::string::gpdType e.g. "Ht"
        param::string::flavor e.g. "dv"
        param::float::x between 0,1 (not the 0 itself)
        param::float::t Negative values (t=0 returns the forward limit)
        """
        nominal = self.xGPD(analysisSet, gpdType, flavor, x, t)
        profFuncParameters = getProfileFunctionParameters(self.__analysis_type, gpdType, analysisSet)(flavor)
        profileFunction = ProfileFunction(profFuncParameters,x)() # last () is intentional, don't remove
        #M = 1 #FIX THIS #FIX THIS
        if "H" == gpdType: ### M,UPDF could be technically outside but would make it a bit complex so let's stick with this
            M = UncertaintyGPD(analysisSet,gpdType, flavor, x,t).uncertainty
            UPDF = [self.__pdfHandler__(flavor,x,self.__UPDF),self.__uncertainPDF__(flavor, x,self.__UGridPDFSet,self.__UGridPDFs)]
            delta = np.sqrt(np.power(np.exp(t*profileFunction) * UPDF[1],2)  +  np.power(UPDF[0]*M , 2))
        if "Ht" == gpdType:
            M = UncertaintyGPD(analysisSet,gpdType, flavor, x,t).uncertainty
            UPDF = [self.__pdfHandler__(flavor,x,self.__PPDF),self.__uncertainPDF__(flavor, x,self.__PGridPDFSet,self.__PGridPDFs)]
            delta = np.sqrt(np.power(np.exp(t*profileFunction) * UPDF[1],2)  +  np.power(UPDF[0]*M , 2))
        if "E" == gpdType:
            delta = UncertaintyGPD(analysisSet,gpdType, flavor, x,t).uncertainty

        return ufloat(nominal,delta)
    
    def xGPDxi(self,analysisSet, gpdType, flavor, x, t, xi):
        """
        param::string::analysisSet e.g. "Set11"
        param::string::gpdType e.g. "Ht"
        param::string::flavor e.g. "dv"
        param::float::x between 0,1 (not the 0 itself)
        param::float::t Negative values (t=0 returns the forward limit)
        param::float::xi between 0,1 (not the 0 itself)
        """
        if xi ==0:
            return self.xGPD(analysisSet, gpdType, flavor, x, t)
        b0 = np.divide(x+xi,1+xi)
        if x <= xi:
            a0 = 1e-5
        else:
            a0 = np.divide(x-xi,1-xi) 
        return quad(self.__xGPDxiIntegrand__, a0, b0, args=(analysisSet, gpdType , flavor,x,t,xi),epsabs=1e-9, limit = 150 )[0]

        
        
    
#################################### Setters
    def list_GPDTypes(self):
        """
        List all directories in the 'src/GPD/data' folder.
        """
        directory = 'src/data/'+self.__analysis_type
        return [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]

    def __get_Q2__(self):
        """
        Return the Q^2 value for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return 4.0
        else:
            return 4.0  # Default value

    def print_analysis_doi(self):
        """
        Print the DOI for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            print("#############################################################################")
            print("Thanks for using our analysis! The corresponding paper is: arXiv:2211.09522v2")
            print("#############################################################################")
            print("\n")

    def __get_analysis_updf__(self):
        """
        Get the unpolarized PDF (UPDF) for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return lhapdf.mkPDF("NNPDF40_nlo_as_01180", 0)

    def __get_analysis_ppdf__(self):
        """
        Get the polarized PDF (PPDF) for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return lhapdf.mkPDF("NNPDFpol11_100", 0)

    def __select_ugrid_pdf_set__(self):
        """
        Select the unpolarized grid PDF set for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return lhapdf.getPDFSet("NNPDF40_nlo_as_01180")

    def __select_pgrid_pdf_set__(self):
        """
        Select the polarized grid PDF set for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return lhapdf.getPDFSet("NNPDFpol11_100")

    def __get_ugrid_pdfs__(self):
        """
        Get the unpolarized grid PDFs for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return self.__UGridPDFSet.mkPDFs()

    def __get_pgrid_pdfs__(self):
        """
        Get the polarized grid PDFs for the given analysis type.
        """
        if self.__analysis_type == "HGAG23":
            return self.__PGridPDFSet.mkPDFs()
################################### End Setters

################################### xGPD Subroutines

    def __pdfHandler__(self, flavor, x,mkPDF):

        sqrtQ = np.sqrt(self.Q2)
        code = self.__flavor_map.get(flavor)
        if isinstance(code, tuple):  # For valence (e.g., "uv", "dv")
            return mkPDF.xfxQ(code[0], x, sqrtQ) - mkPDF.xfxQ(code[1], x, sqrtQ)
        elif code is not None:  # For other flavors (e.g., "u", "ubar")
            return mkPDF.xfxQ(code, x, sqrtQ)
        else:
            raise ValueError(f"Unknown flavor: {flavor}")
        
    def __uncertainPDF__(self, flavor, x,pset,mkpdf):
        sqrtQ = np.sqrt(self.Q2)
        code = self.__flavor_map.get(flavor)
        xfAll = [0.0] * pset.size
        # Fill xfAll based on flavor 
        if isinstance(code, tuple):  # For valence (e.g., "uv", "dv")
            for imem in range(pset.size):
                xfAll[imem] = mkpdf[imem].xfxQ(code[0], x, sqrtQ) - mkpdf[imem].xfxQ(code[1], x, sqrtQ)
        elif code is not None:  # For other flavors (e.g., "u", "ubar")
            for imem in range(pset.size):
                xfAll[imem] = mkpdf[imem].xfxQ(code, x, sqrtQ)
        else:
            raise ValueError(f"Unknown flavor: {flavor}")

        # Compute uncertainty
        Uncf = pset.uncertainty(xfAll, cl=pset.errorConfLevel)
        w = (Uncf.errplus + Uncf.errminus) / 2  # Average uncertainty
        s = Uncf.scale
        return w * (1 / s)



    def __pdfEHandler__(self,analysisSet, flavor, x):
        def integrand(x, alpha_q, beta_q, gamma_q):
            return (np.power(x, -alpha_q) * np.power(1-x, beta_q) * (1 + gamma_q * np.sqrt(x)))
        k = {
        "uv": 1.67,
        "dv": -2.03
        }
        parameterList = getProfileFunctionParameters(self.__analysis_type, "E", analysisSet)(flavor) 
        N = np.divide(1, quad(integrand, 0, 1, args=(parameterList[3], parameterList[4], parameterList[5]))[0])
        return x * k.get(flavor) * N * np.power(x,-parameterList[3]) * np.power(1-x,parameterList[4]) * (1+ parameterList[5] * np.sqrt(x))
    

    def __xGPDxiIntegrand__( self, b, analysisSet, gpdType , flavor , x, t , xi ):
        #Hv = MMGPD.xGPD(InitilizerArgs, Set, GPDType , Flavour, b, t) / b
        Hv = self.xGPD(analysisSet, gpdType, flavor, x, t)
        sd = np.divide(Hv, np.power(1-b,3))
        return np.divide(3,4)*sd*(np.power(1-b,2)-np.power(x-b,2)/np.power(xi,2))/xi

################################### End of xGPD Subroutines
