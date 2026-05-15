#from typing_extensions import ParamSpecArgs
import sympy as sp
import numpy as np

# : a changer en cas de problem (formule not just peut être)


class Ply:
    def __init__(self, thickness, fiberMaterial, fiberVolume, fiberOrientation, matrixMaterial, e2rule="Rule Of Mixture", g12rule="Rule Of Mixture"):
        self.t = thickness
        self.fM = fiberMaterial
        self.Vf = fiberVolume
        self.tetha = fiberOrientation  # in degree
        self.mM = matrixMaterial
        self.Vm = 1 - fiberVolume  # as Vf + Vm = 1 for each ply
        # Elastic constants
        # longitidunal stiffness(Young's modulus of the ply in the 1-direction)
        self.e1 = round(fiberMaterial.e1 * self.Vf +
                        matrixMaterial.e * self.Vm, 4) # CHECKED
        self.v12 = round(fiberMaterial.v12 * self.Vf + matrixMaterial.v *
                         self.Vm, 4)  # Poisson's ratio v12  ("nu" 1-2) of the Ply #CHECKED
        # CHECKED (the default method,modified ROM (1&2))
        self.e2 = self.e2(e2rule)
        self.g12 = self.g12(g12rule)  # CHECKED (for all 3 methods)
        """self.e1 = 155  # GPa
        self.e2 = 12.10  # GPa
        self.v12 = 0.248
        self.g12 = 4.40  # GPa"""

    # transverse stiffness (Young's modulus in the 2-direction)

    def e2(self, method="Rule Of Mixture"):  # CHECKED (ALL METHODS)
        e2 = 0
        if method == "Rule Of Mixture":
            em_prime = self.mM.e/(1-self.mM.v)
            e2 = round(1 / (self.Vf/self.fM.e2 +
                            self.Vm/em_prime), 4)  # CHECKED
        elif method == "Halpin-Tsai":
            # Rigorously zeta is determinated through experimental data, but here we use the value of 2 (Katchy,2008)
            zeta = 2
            etha = (self.fM.e1-self.mM.e)/(self.fM.e1+zeta*self.mM.e)
            e2 = self.mM.e*((1+zeta*etha*self.Vf)/(1-etha*self.Vf))
        elif method == "Voyiadjis-Kattan1":
            etha = 0.5  # it is usually taken between 0.4 and 0.6
            e2 = (self.Vf+etha*self.Vm) / \
                (self.Vf/self.fM.e2 + etha*self.Vm/self.mM.e)  # CHECKED
        elif method == "Voyiadjis-Kattan2":
            e1f = self.fM.e1
            e2f = self.fM.e2
            em = self.mM.e
            vf = self.Vf
            vm = self.Vm
            nu12f = self.fM.v12
            num = self.mM.v
            ethaf = (e1f*vf + ((1-nu12f**2)*em + num*nu12f*e1f)*vm) / \
                (e1f*vf+em*vm)
            etham = (((1-num**2)*e1f - (1-num*nu12f)*em)
                     * vf + em*vm)/(e1f*vf + em*vm)
            e2 = e2f*em/(ethaf*vf*em + etham*vm*e2f)  # CHECKED
        else:
            raise("Unkown method for  ply E2")

        return e2

    def g12(self, method="Rule Of Mixture"):  # CHECKED (ALL METHODS)
        g12 = 0
        g12f = self.fM.g12
        gm = self.mM.g
        vf = self.Vf
        vm = self.Vm
        if method == "Rule Of Mixture":
            g12 = round(1/(self.Vf/self.fM.g12 + self.Vm /
                           self.mM.g), 4)  # CHECKED
        elif method == "Halpin-Tsai":
            ksi = 1  # best agreement with experimental results
            etha = (self.fM.g12 - self.mM.g)/(self.fM.g12 + ksi*self.mM.g)
            g12 = round(self.mM.g*((1+ksi*etha*self.Vf)/(1-etha*self.Vf)), 4)
        elif method == "modified-ruleOfMixture":
            # this value of ethaPrime leads to good correlation with elasticity solution (refxx)
            ethaPrime = 0.6
            g12 = 1/((vf/g12f + ethaPrime*vm/gm)/(vf+ethaPrime*vm))
        elif method == "self-consistent field":
            g12 = self.gm*(((1+self.Vf)*self.fM.g12+self.Vm*self.mM.g) /
                           (self.Vm*self.fM.g12 + (1+self.Vf)*self.mM.g))
        else:
            raise("Unknow method for G12")
        return g12

    def sMatrix(self):  # REDUCED COMPLIANCE MATRIX CHECKED
        s = np.array([[1/self.e1, -self.v12/self.e1, 0],
                      [-self.v12/self.e1, 1/self.e2, 0],
                      [0, 0, 1/self.g12]])
        return s

    def qMatrix(self):  # REDUCED STIFFNESS MATRIX CHECKED
        return np.linalg.inv(self.sMatrix())

    def tMatrix(self):
        tta = self.tetha*np.pi/180  # fiber orientation in radian
        cos2 = np.cos(tta)**2  # COS^2(thetha)
        sin2 = np.sin(tta)**2  # SIN^2(tetha)
        cos = np.cos(tta)  # COS(tetha)
        sin = np.sin(tta)  # SIN(tetha)
        t = np.array([[cos2, sin2, 2*cos*sin],
                      [sin2, cos2, -2*cos*sin],
                      [-sin*cos, sin*cos, cos2-sin2]])
        return t

    def sbarMatrix(self):  # CHECKED
        t = self.tMatrix()
        s = self.sMatrix()
        s[2][2] = 0.5*s[2][2]
        t_inv = np.linalg.inv(t)
        sbar = np.matmul(np.matmul(t_inv, s), t)
        for i in range(3):
            sbar[2][i] = 2*sbar[2][i]
        """for i in range(len(sbar)):
            for j in range(len(sbar[i])):
                sbar[i][j] = round(sbar[i][j], 4)"""
        return sbar

    def qbarMatrix(self):  # CHECKED
        sbar = self.sbarMatrix()
        qbar = np.linalg.inv(sbar)
        """for i in range(len(qbar)):
            for j in range(len(qbar[i])):
                qbar[i][j] = round(qbar[i][j], 4)"""
        return qbar
