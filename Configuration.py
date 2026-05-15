from Ply import Ply
from Laminate import Laminate
import pandas as pd
import numpy as np
import sympy as sp
from materials import Fiber, Matrix
import os


class Configuration:
    def __init__(self, excel_laminate_data):
        self.dfname = excel_laminate_data
        self.laminate_df = pd.read_excel(
            self.dfname, index_col=0)
        self.orientations = np.array(
            self.laminate_df["Fibers\nOrientation"].tolist())

        self.fibers = np.array(
            self.laminate_df["Fiber"].tolist())

        self.fiberVolumes = np.array(
            self.laminate_df["Fiber\n Volume"].tolist())

        self.matrices = np.array(self.laminate_df['Matrix'].tolist())

        self.thicknesses = np.array(
            self.laminate_df["Ply thickness (mm)"].tolist())
        self.numberOfPlies = len(self.orientations)

    def plies(self):
        plies = [Ply(self.thicknesses[i], Fiber(self.fibers[i]), self.fiberVolumes[i],
                     self.orientations[i], Matrix(self.matrices[i])) for i in range(self.numberOfPlies)]
        return plies

    def forces(self):  # in Newton
        dfForces = pd.read_excel(self.dfname, "Forces and Moments")
        Nx = dfForces["Nx (N)"].tolist()[0]
        Ny = dfForces["Ny (N)"].tolist()[0]
        Nxy = dfForces["Nxy (N)"].tolist()[0]
        return Nx, Ny, Nxy

    def moments(self):  # in Newton/meter
        dfMoments = pd.read_excel(self.dfname, "Forces and Moments")
        Mx = dfMoments["Mx (N.m)"].tolist()[0]
        My = dfMoments["My (N.m)"].tolist()[0]
        Mxy = dfMoments["Mxy (N.m)"].tolist()[0]
        return Mx, My, Mxy

    # in Newton/meter (assuming uniforme distribution of Force per unit of length)
    def forcesPerUnitLength(self):
        dfForces = pd.read_excel(self.dfname, "Forces and Moments")
        dflength = pd.read_excel(self.dfname, "Laminate conf")
        Lx = dflength["Lx (m)"].tolist()[0]
        Ly = dflength["Ly (m)"].tolist()[0]
        Nx = dfForces["Nx (N)"].tolist()[0]/Ly
        Ny = dfForces["Ny (N)"].tolist()[0]/Lx
        # With the line 50 assumption, Nxy not equal Nyx but Nxy/Ly = Nyx/Lx (for equilibrium)
        Nxy = dfForces["Nxy (N)"].tolist()[0]/Ly

        return Nx, Ny, Nxy

    def momentsPerUnitLenght(self):
        dfMoments = pd.read_excel(self.dfname, "Forces and Moments")
        dflength = pd.read_excel(self.dfname, "Laminate conf")
        Lx = dflength["Lx (m)"].tolist()[0]
        Ly = dflength["Ly (m)"].tolist()[0]
        Mx = dfMoments["Mx (N.m)"].tolist()[0]/Ly
        My = dfMoments["My (N.m)"].tolist()[0]/Lx
        # see explications in line 58 :"With the line 50 assumption..."
        Mxy = dfMoments["Mxy (N.m)"].tolist()[0]/Ly
        return Mx, My, Mxy
