import numpy as np
import os
#from numpy.core.defchararray import lower

from numpy.lib.function_base import select
from PlyStrength import PlyStrength

"""The commentary #CHECKED means that function was verified by testing with numerical values from book
 # CHECKED for symmetric means the function is checked only for symmetric laminate """


class Laminate:
    def __init__(self, plies):
        self.plies = plies
        self.qbarMatrices = np.array(
            [self.plies[i].qbarMatrix() for i in range(len(plies))])
        self.ts = np.array([plies[i].t for i in range(len(plies))])
        self.thickness = sum(self.ts)
        print("thickness: ", self.thickness)
        # zk is the list containing the different height of the interfaces
        zk = [-self.thickness/2]
        # filling zk
        for i in range(len(plies)):
            zk.append(round(zk[i]+self.ts[i], 4))
        self.zk = zk  # du dernier au premier
        print("les zk")
        print(zk)

    def AMatrix(self):  # CHECKED
        A = list()
        for i in range(3):
            line_A = [sum([self.qbarMatrices[k][i][j] * self.ts[k]
                           for k in range(len(self.qbarMatrices))]) for j in range(3)]
            A.append(line_A)
        A = np.array(A)
        return A

    def BMatrix(self):  # CHECKED
        B = list()
        for i in range(3):
            line_B = [0.5*sum([self.qbarMatrices[k][i][j] * (self.zk[k+1]**2 - self.zk[k]**2)
                               for k in range(len(self.qbarMatrices))])for j in range(3)]
            B.append(line_B)
        B = np.array(B)
        return B

    def DMatrix(self):  # CHECKED
        D = list()
        for i in range(3):
            line_d = [1/3 * sum([self.qbarMatrices[k][i][j] * (self.zk[k+1]**3 - self.zk[k]**3)
                                 for k in range(len(self.qbarMatrices))])for j in range(3)]
            D.append(line_d)
        D = np.array(D)
        return D

    def ABDMatrix(self):  # checked since abcd is checked
        """return ABD Matrix by concatenating A, B and D Matrices"""
        A = self.AMatrix() * 10**6  # in Newton/meter
        B = self.BMatrix()*10**3  # in Newton
        D = list(self.DMatrix())  # in Newton.meter
        AB = np.concatenate((A, B), axis=1)
        BD = np.concatenate((B, D), axis=1)
        ABD = np.concatenate((AB, BD), axis=0)
        return ABD

    def abcdMatrix(self):  # CHECKED with examples in stress and strain/ Forces and moments
        ABD = self.ABDMatrix()
        return np.linalg.inv(ABD)

    def aMatrix(self):
        """this function extract a_Matrix from abcd one"""
        abcd = self.abcdMatrix()
        ab, cd = np.vsplit(abcd, 2)
        a, b = np.hsplit(ab, 2)
        c, d = np.hsplit(cd, 2)
        """Direct computation method:
        A = self.AMatrix()
        B = self.BMatrix()
        D = self.DMatrix()
        Ainv = np.linalg.inv(A)
        B_start = np.matmul(-Ainv, B)
        C_start = np.matmul(B, Ainv)
        Dstart = D - np.matmul(np.matmul(B, Ainv), B)
        Dstartinv = np.linalg.inv(Dstart)
        a = Ainv - np.matmul(np.matmul(B_start, Dstartinv), C_start)"""
        return a

    def bMatrix(self):
        """this function extract b_Matrix from abcd one"""
        abcd = self.abcdMatrix()
        ab, cd = np.vsplit(abcd, 2)
        a, b = np.hsplit(ab, 2)
        c, d = np.hsplit(cd, 2)
        """Direct computation method:
        A = self.AMatrix()
        B = self.BMatrix()
        D = self.DMatrix()
        Ainv = np.linalg.inv(A)
        B_start = np.matmul(-Ainv, B)
        C_start = np.matmul(B, Ainv)
        Dstart = D - np.matmul(np.matmul(B, Ainv), B)
        Dstartinv = np.linalg.inv(Dstart)
        b = np.matmul(B_start, Dstartinv)
        """
        return b

    def cMatrix(self):
        """this function extract c_Matrix from abcd one"""
        abcd = self.abcdMatrix()
        ab, cd = np.vsplit(abcd, 2)
        a, b = np.hsplit(ab, 2)
        c, d = np.hsplit(cd, 2)
        """Direct computation method:
        A = self.AMatrix()
        B = self.BMatrix()
        D = self.DMatrix()
        Ainv = np.linalg.inv(A)
        B_start = np.matmul(-Ainv, B)
        C_start = np.matmul(B, Ainv)
        Dstart = D - np.matmul(np.matmul(B, Ainv), B)
        Dstartinv = np.linalg.inv(Dstart)
        c = np.matmul(-Dstartinv, C_start)"""
        return c

    def dMatrix(self):
        """this function extract d_Matrix from abcd one"""
        abcd = self.abcdMatrix()
        ab, cd = np.vsplit(abcd, 2)
        a, b = np.hsplit(ab, 2)
        c, d = np.hsplit(cd, 2)
        """Direct computation method:
        A = self.AMatrix()
        B = self.BMatrix()
        D = self.DMatrix()
        Ainv = np.linalg.inv(A)
        B_start = np.matmul(-Ainv, B)
        C_start = np.matmul(B, Ainv)
        Dstart = D - np.matmul(np.matmul(B, Ainv), B)
        d = np.linalg.inv(Dstart)"""
        return d

    def elastic_constants(self):  # CHECKED FOR SYMMETRIC LAMINATE
        """This function return the elastic constants of the laminate (Ex, Ey, Gxy, vxy, vyx)"""
        a = 10**9*self.aMatrix()
        Ex = 1/((a[0][0])*self.thickness*0.001)
        Ey = 1/((a[1][1])*self.thickness*0.001)
        Gxy = 1/((a[2][2])*self.thickness*0.001)
        vxy = -a[1][0]/a[0][0]
        vyx = -a[0][1]/a[1][1]
        return round(Ex, 4), round(Ey, 4), round(Gxy, 4), round(vxy, 4), round(vyx, 4)

    def stressesInLamina_interfaces(self, midplane_strain):
        """This function compute sigma_x,sigma_y,Tau_xy in the upper,lower and middle z coordinates
        of each ply and print them"""

        # laminate different interfaces coordinate (lower and upper z coordinate of each lamina)
        zk = self.zk

        # Computation for each lamina (or ply)
        for i in range(len(self.plies)):
            qbarmatrice = self.qbarMatrices[i]  # GPa

            # lamina upper bound z coordinate
            zu = zk[i]*10**-3  # meter
            # then, lamina upper bound stresses
            upperStresses = np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zu*np.matmul(qbarmatrice, midplane_strain[3:])

            # lamina middle z coordinate
            zm = ((zk[i]+zk[i+1])/2)*10**-3  # meter
            # then, lamina middle bound stresses
            middleStresses = np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zm*np.matmul(qbarmatrice, midplane_strain[3:])

            # lamina lower bound z coordinate
            zl = zk[i+1]*10**-3  # meter
            # then, lamina lower bound stresses
            lowerStresses = np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zl*np.matmul(qbarmatrice, midplane_strain[3:])

            # Now printing  results :-)
            print("ply ", i+1, ": ", self.plies[i].tetha, "°:\n")
            print("zu: ", zu*1000, "mm: sigma_x: ", round(upperStresses[0][0], 4), " GPa\t", "sigma_y: ",
                  round(upperStresses[1][0], 4), " GPa\t", "taux_xy: ", round(upperStresses[2][0], 4), " GPa")
            print("zm: ", zm*1000, "mm: sigma_x: ", round(middleStresses[0][0], 4), " GPa\t", "sigma_y: ",
                  round(middleStresses[1][0], 4), " GPa\t", "taux_xy: ", round(middleStresses[2][0], 4), " GPa")
            print("zl: ", zl*1000, "mm: sigma_x: ", round(lowerStresses[0][0], 4), " GPa\t", "sigma_y: ",
                  round(lowerStresses[1][0], 4), " GPa\t", "taux_xy: ", round(lowerStresses[2][0], 4), " GPa", end="\n\n")

    def principallAxesStressesInLamina_interfaces(self, midplane_strain):
        """This function compute sigma1,sigman2,Tau12 in the upper,lower and middle z coordinate
        for each ply and print them"""

        # laminate different interfaces coordinate (lower and upper z coordinate of each lamina)
        zk = self.zk

        # Computation for each lamina (or ply)
        for i in range(len(self.plies)):
            qbarmatrice = self.qbarMatrices[i]
            tMatrix = self.plies[i].tMatrix()

            # lamina upper bound z coordinate
            zu = zk[i]*10**-3  # meter
            # then, lamina upper bound stresses
            upperStresses = np.matmul(tMatrix, np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zu*np.matmul(qbarmatrice, midplane_strain[3:]))

            # lamina middle z coordinate
            zm = (zk[i]+zk[i+1])/2*10**-3  # meter
            # then, lamina middle bound stresses
            middleStresses = np.matmul(tMatrix, np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zm*np.matmul(qbarmatrice, midplane_strain[3:]))

            # lamina lower bound z coordinate
            zl = zk[i+1]*10**-3  # meter
            # then, lamina lower bound stresses
            lowerStresses = np.matmul(tMatrix, np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zl*np.matmul(qbarmatrice, midplane_strain[3:]))

            # printing results
            print("ply ", i+1, ": ", self.plies[i].tetha, "°:\n")
            print("zu: ", zu*1000, "mm: sigma_1: ", round(upperStresses[0][0], 6), " GPa\t", "sigma_2: ",
                  round(upperStresses[1][0], 6), " GPa\t", "taux_12: ", round(upperStresses[2][0], 6), " GPa")
            print("zm: ", zm*1000, "mm: sigma_1: ", round(middleStresses[0][0], 6), " GPa\t", "sigma_2: ",
                  round(middleStresses[1][0], 6), " GPa\t", "taux_12: ", round(middleStresses[2][0], 6), " GPa")
            print("zl: ", zl*1000, "mm: sigma_1: ", round(lowerStresses[0][0], 6), " GPa\t", "sigma_2: ",
                  round(lowerStresses[1][0], 6), " GPa\t", "taux_12: ", round(lowerStresses[2][0], 6), " GPa", end="\n\n")

    def tsaiHill(self, midplane_strain):
        # evaluated for upper and lower  bounds of each lamina (the stresses are maximum)
        zk = self.zk
        for i in range(len(self.plies)):
            qbarmatrice = self.qbarMatrices[i]
            tMatrix = self.plies[i].tMatrix()
            zu = zk[i]*10**-3  # meter  # lamina upper bound z coordinate
            upperStresses = np.matmul(tMatrix, np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zu*np.matmul(qbarmatrice, midplane_strain[3:]))
            zl = zk[i+1]*10**-3  # meter  # lamina lower bound z coordinate
            lowerStresses = np.matmul(tMatrix, np.matmul(
                qbarmatrice, midplane_strain[0:3]) + zl*np.matmul(qbarmatrice, midplane_strain[3:]))
            print("ply ", i+1, ": ", self.plies[i].tetha, "°:\n")
            # considered ply strength: "c" for compression "t" for tension and "s" for shear
            strength = PlyStrength(self.plies[i])
            F1t = strength.tensileStrenght1()*10**-3  # GPa
            F2t = strength.tensileStrength2()*10**-3  # GPa
            F1c = strength.compressiveStrength1()*10**-3  # GPa
            F2c = strength.compressiveStrength2()*10**-3  # GPa
            Fs = strength.shearStrength()*10**-3  # GPa
            alfasquared1 = 0
            alfasquared2 = 0
            alfasquared3 = 0
            alfasquared4 = upperStresses[2][0]**2/Fs**2
            # upper stresses
            if upperStresses[0][0] >= 0:  # tensile stress
                alfasquared1 = upperStresses[0][0]**2/F1t**2
                if upperStresses[1][0] >= 0:  # if sigma2>0
                    alfasquared2 = upperStresses[0][0] * \
                        upperStresses[1][0]/F1t**2
                    alfasquared3 = upperStresses[1][0]**2/F2t**2
                else:
                    alfasquared2 = upperStresses[0][0] * \
                        upperStresses[1][0]/(F1t*F1c)
                    alfasquared3 = upperStresses[1][0]**2/F2c**2

            else:
                alfasquared1 = upperStresses[0][0]**2/F1c**2
                if upperStresses[1][0] >= 0:  # if sigma2>0
                    alfasquared2 = upperStresses[0][0] * \
                        upperStresses[1][0]/F1c*F1t
                    alfasquared3 = upperStresses[1][0]**2/F2t**2
                else:
                    alfasquared2 = upperStresses[0][0] * \
                        upperStresses[1][0]/(F1c**2)
                    alfasquared3 = upperStresses[1][0]**2/F2c**2
            alfasquared = alfasquared1-alfasquared2+alfasquared3+alfasquared4
            if alfasquared <= 1:
                print("lamina ", i+1,
                      " (", self.plies[i].tetha, "°):\nupper bound OK")
            else:
                print("failure in lamina ", i+1,
                      " (", self.plies[i].tetha, "°): upper bound\n")
            # lower stresses
            if lowerStresses[0][0] >= 0:  # tensile stress
                alfasquared1 = lowerStresses[0][0]**2/F1t**2
                if lowerStresses[1][0] >= 0:  # if sigma2>0
                    alfasquared2 = lowerStresses[0][0] * \
                        lowerStresses[1][0]/F1t**2
                    alfasquared3 = lowerStresses[1][0]**2/F2t**2
                else:
                    alfasquared2 = lowerStresses[0][0] * \
                        lowerStresses[1][0]/(F1t*F1c)
                    alfasquared3 = lowerStresses[1][0]**2/F2c**2

            else:
                alfasquared1 = lowerStresses[0][0]**2/F1c**2
                if lowerStresses[1][0] >= 0:  # if sigma2>0
                    alfasquared2 = lowerStresses[0][0] * \
                        lowerStresses[1][0]/F1c*F1t
                    alfasquared3 = lowerStresses[1][0]**2/F2t**2
                else:
                    alfasquared2 = lowerStresses[0][0] * \
                        lowerStresses[1][0]/(F1c**2)
                    alfasquared3 = lowerStresses[1][0]**2/F2c**2
            alfasquared = alfasquared1-alfasquared2+alfasquared3+alfasquared4
            if alfasquared <= 1:
                print(print("lamina ", i+1,
                            " (", self.plies[i].tetha, "°):\nlower bound OK"))
            else:
                print("failure lamina ", i+1,
                      " (", self.plies[i].tetha, " °): lower bound\n")

    def axialStrength(self, stressType="1_tensile,2_tensile"):
        """this function computes the axial strength (sigma_x max) for each ply and print them"""

        for i in range(len(self.plies)):
            # concerned ply orientation in radian
            tetha_rad = self.plies[i].tetha*np.pi/180
            # cos^4,cos^2,sin^4,sin^2
            m4, m2, n4, n2 = np.cos(
                tetha_rad)**4, np.cos(tetha_rad)**2, np.sin(tetha_rad)**4, np.sin(tetha_rad)**2
            strength = PlyStrength(self.plies[i])
            F1t = strength.tensileStrenght1()*10**-3  # GPa
            F2t = strength.tensileStrength2()*10**-3  # GPa
            F1c = strength.compressiveStrength1()*10**-3  # GPa
            F2c = strength.compressiveStrength2()*10**-3  # GPa
            Fs = strength.shearStrength()*10**-3  # GPa
            Fx = 0
            if stressType == "1_tensile,2_tensile":
                Fx = np.sqrt(1/(m4/F1t**2 + n4/F2t**2 +
                                (1/Fs**2 - 1/F1t**2)*m2*n2))
            elif stressType == "1_tensile,2_compressive":
                Fx = np.sqrt(1/(m4/F1t**2 + n4/F2c**2 +
                                (1/Fs**2 - 1/F1t**2)*m2*n2))
            elif stressType == "1_compressive,2_tensile":
                Fx = np.sqrt(1/(m4/F1c**2 + n4/F2t**2 +
                                (1/Fs**2 - 1/F1c**2)*m2*n2))
            elif stressType == "1_compressive,2_compressive":
                Fx = np.sqrt(1/(m4/F1c**2 + n4/F2c**2 +
                                (1/Fs**2 - 1/F1c**2)*m2*n2))
            else:
                raise("Unknown sressType for axial Strength (sigma_x Max) function")
            print("ply ", i+1, ": ", self.plies[i].tetha, "°:\nSigmax_Max:")
            print(round(Fx*1000, 6), "MPa\n\n")
