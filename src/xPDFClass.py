import numpy as np
import lhapdf
from uncertainties import ufloat

class xPDF:
    def __init__(self, pdfName):
        self.cen = lhapdf.mkPDF(pdfName, 0)  # Central PDF member
        self.pset = lhapdf.getPDFSet(pdfName)  # PDF set
        self.pdfs = self.pset.mkPDFs()  # All PDF members
        self.xfAll = [0.0] * self.pset.size  # Initialize xfAll array

        # Map flavors to LHAPDF codes
        self.flavor_map = {
            "u": 2, "ubar": -2, "uv": (2, -2),
            "d": 1, "dbar": -1, "dv": (1, -1),
            "s": 3, "sbar": -3, "sv": (3, -3),
            "c": 4, "cbar": -4, "cv": (4, -4),
            "b": 5, "bbar": -5, "bv": (5, -5),
            "t": 6, "tbar": -6, "tv": (6, -6),
            "g": 21
        }

    def xPDFwUncertinty(self, flavor, x, Q2):
        return ufloat(self.xPDFCentVal(flavor, x, Q2), self.__uncertainty__(flavor, x, Q2))

    def xPDFCentVal(self, flavor, x, Q2):
        sqrtQ = np.sqrt(Q2)
        code = self.flavor_map.get(flavor)

        if isinstance(code, tuple):  # For valence (e.g., "uv", "dv")
            return self.cen.xfxQ(code[0], x, sqrtQ) - self.cen.xfxQ(code[1], x, sqrtQ)
        elif code is not None:  # For other flavors (e.g., "u", "ubar")
            return self.cen.xfxQ(code, x, sqrtQ)
        else:
            raise ValueError(f"Unknown flavor: {flavor}")

    def __uncertainty__(self, flavor, x, Q2):
        sqrtQ = np.sqrt(Q2)
        code = self.flavor_map.get(flavor)

        # Fill xfAll based on flavor
        if isinstance(code, tuple):  # For valence (e.g., "uv", "dv")
            for imem in range(self.pset.size):
                self.xfAll[imem] = self.pdfs[imem].xfxQ(code[0], x, sqrtQ) - self.pdfs[imem].xfxQ(code[1], x, sqrtQ)
        elif code is not None:  # For other flavors (e.g., "u", "ubar")
            for imem in range(self.pset.size):
                self.xfAll[imem] = self.pdfs[imem].xfxQ(code, x, sqrtQ)
        else:
            raise ValueError(f"Unknown flavor: {flavor}")

        # Compute uncertainty
        Uncf = self.pset.uncertainty(self.xfAll, cl=self.pset.errorConfLevel)
        w = (Uncf.errplus + Uncf.errminus) / 2  # Average uncertainty
        s = Uncf.scale
        return w * (1 / s)
