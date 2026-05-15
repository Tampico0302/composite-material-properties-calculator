import pandas as pd
import os


class Fiber:
    def __init__(self, fiberName):
        # fiberName is of str type
        df_Fiber = pd.read_excel(
            "settings.xlsx", sheet_name="Fibers", index_col=0)
        self.fiberName = fiberName

        self.e1 = df_Fiber["El(GPa)"][fiberName]
        self.e2 = df_Fiber["Et(GPa)"][fiberName]
        self.v12 = df_Fiber["vlt"][fiberName]
        self.g12 = df_Fiber["Glt-(GPa)"][fiberName]
        # Tensile strength (fiber direction)
        self.Ft = df_Fiber["Flt-(Mpa)"][fiberName]


class Matrix:
    def __init__(self, matrixName):
        df_Matrix = pd.read_excel(
            "settings.xlsx", sheet_name="Matrix", index_col=0)
        self.matrixName = matrixName
        self.e = df_Matrix["E (GPa)"][matrixName]
        self.g = df_Matrix["G (GPa)"][matrixName]
        self.v = df_Matrix["v"][matrixName]
        self.Ft = df_Matrix["Ft (MPa)"][matrixName]  # Tensile strength
        self.Fc = df_Matrix["Fc (MPa)"][matrixName]  # compressive strength
        self.Fs = df_Matrix["Fms (MPa)"][matrixName]  # shear strength
