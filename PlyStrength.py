import numpy as np


class PlyStrength:
    def __init__(self, ply):
        # ultimate fiber tensile strain (brittle fiber)
        self.uftstr = ply.fM.Ft/ply.fM.e1
        self.umtstr = ply.mM.Ft
        self.fFt = ply.fM.Ft  # fiber tensile Strength
        self.mFt = ply.mM.Ft  # matrix tensile Strength
        self.mFc = ply.mM.Fc  # matrix comrpressive Strength
        self.mFs = ply.mM.Fs  # matrix shear Strength
        self.em = ply.mM.e
        self.ef1 = ply.fM.e1
        self.Vf = ply.Vf
        self.Vm = ply.Vm
        self.gm = ply.mM.g
        g12 = ply.g12
        self.e1 = ply.e1
        self.nu_m = ply.mM.v
        self.g12f = ply.fM.g12

    def tensileStrenght1(self):
        Ft1 = 0  # initialisation of tensile strength
        if self.uftstr < self.umtstr:
            VminFiber = (self.mFt - self.em * self.uftstr) / \
                (self.mFt + self.fFt - self.em*self.uftstr)
            """VcritFiber = (self.mFt - self.em * self.uftstr) / \
                (self.fFt - self.em*self.uftstr)"""
            if self.Vf < VminFiber:
                Ft1 = self.mFt*self.Vm
            else:
                Ft1 = self.fFt*self.Vf + self.em*self.uftstr*self.Vm
        else:
            Ft1 = self.mFt*(self.Vf*self.ef1/self.em + self.Vm)
        return Ft1

    def compressiveStrength1(self, method="rosen shear mode"):
        Fc1 = 0
        if method == "rosen shear mode":
            Fc1 = self.gm/self.Vm
        elif method == "rosen extensional mode":
            Fc1 = 2*self.Vf*(self.em*self.ef1*self.Vf/(3*self.Vm))**0.5
        elif method == 'Lo-chim model':  # good correlation for Eglass epoxy
            Fc1 = self.g12/(1.5+12*(6/np.pi)**2 * (self.g12/self.e1))
        else:
            raise("Unkown compressive strenght prediction method")
        return Fc1

    def tensileStrength2(self, method="maximum tensile stress criterion"):
        """residual stresses are not considered"""
        Ft2 = 0
        ksigma = (1-self.Vf*(1-self.em/self.ef1)) / \
            (1-(4*self.Vf/np.pi)**0.5 * (1-self.em/self.ef1))
        if method == "maximum tensile stress criterion":
            Ft2 = 1/ksigma * (self.mFt)
        elif method == "maximum tensile stain criterion":
            Ft2 = (1-self.nu_m)*self.mFt/(ksigma*(1+self.nu_m)*(1-2*self.nu_m))
        else:
            raise("Unknown method for transverse tensile strength")
        return Ft2

    def compressiveStrength2(self):
        ksigma = (1-self.Vf*(1-self.em/self.ef1)) / \
            (1-(4*self.Vf/np.pi)**0.5 * (1-self.em/self.ef1))
        return self.mFc/ksigma

    def shearStrength(self):
        ktau = (1-self.Vf*(1-self.gm/self.g12f)) / \
            (1-(4*self.Vf/np.pi)**0.5 * (1-self.gm/self.g12f))
        Fs = self.mFs/ktau
        return Fs
