#from numpy.lib.function_base import _place_dispatcher
from Configuration import Configuration
from Laminate import Laminate
import numpy as np
import pandas as pd
import os

laminateConfiguration = Configuration("settings.xlsx")
plies = laminateConfiguration.plies()
laminate = Laminate(plies)
zk = laminate.zk  # interfaces z coordinates
abcdMatrix = laminate.abcdMatrix()
Nx, Ny, Nxy = laminateConfiguration.forcesPerUnitLength()
Mx, My, Mxy = laminateConfiguration.momentsPerUnitLenght()
forces_moments = np.array([[Nx], [Ny], [Nxy], [Mx], [My], [Mxy]])
laminateMidplanestrains = np.matmul(abcdMatrix, forces_moments)
# laminateMidplanestrains = [[0.001], [0], [0], [0], [0], [0]]


def firstScreen():
    print("input forces (Newton):")
    print("Nx, Ny,Nxy")
    print(laminateConfiguration.forces())
    print("input moments (Newton/meter)")
    print("Mx, My, Mxy")
    print(laminateConfiguration.moments())
    print("-----------------------------------------------------------")
    print("forces per unit of length (N/m")
    print("Nx, Ny, Nxy")
    print(laminateConfiguration.forcesPerUnitLength())
    print("Moments per unit of length (N/m^2)")
    print("Mx,My,Mxy")
    print(laminateConfiguration.momentsPerUnitLenght())
    response = input("Confirmation? (y/n)")
    return response


def secondScreen():
    os.system("cls")
    print("*************input laminate***********")
    print('Number of plies:', laminateConfiguration.numberOfPlies)
    print("plies properties:")
    print("N°  /Fiber\t    /Matrix\t        /tetha  /Vf  /thickness (mm)")
    row = np.array([[i+1, "  /", plies[i].fM.fiberName, "  /", plies[i].mM.matrixName, "  /",
                     plies[i].tetha, "  /", plies[i].Vf, "  /", plies[i].t] for i in range(len(plies))])
    for i in range(len(row)):
        for j in range(len(row[i])):
            print(row[i][j], end=" ")
        print(" ")
    response = input("Confirmation? (y/n)")
    return response


def thirdScreen():
    os.system("cls")
    print("1.)  Engineering elatic constants")
    print("Engineering constants")
    print("Ex(unit), Ey(unit), Gxy(unit), vxy, vyx")
    print(laminate.elastic_constants())
    print("\n***********************************\n")
    print("2. )Stresses (in each laminate) :")
    # reference surface point, strains: epsilonx_0,epsilony_0,gammaxy_0,kx_0,ky_0,kxy_0
    print("\n********midplane strains*******\n")
    print("epsilonx_0,  epsilony_0,  gammaxy_0,  kx_0,  ky_0,  kxy_0:", end="\n\n")
    print([round(laminateMidplanestrains[i][0], 6)
           for i in range(len(laminateMidplanestrains))])
    print("\n\nfor each lamina, stresses in interfaces (upper:zu and lower:zl bound) and middle plane:zm")
    laminate.stressesInLamina_interfaces(laminateMidplanestrains)
    print("\n\nStress in the principal axes for each lamina")
    laminate.principallAxesStressesInLamina_interfaces(laminateMidplanestrains)
    print("\n***********************************\n")
    print("3.)Tsai Hill Expressions")
    print("alpha computation")
    laminate.tsaiHill(laminateMidplanestrains)
    print("\n***********************************\n")
    print("4. ) Maximal allowed stress for each lamina (lamina Strength)")
    for i in range(len(plies)):
        print("ply ", i+1)
        laminate.axialStrength("1_tensile,2_tensile")


def forcingGoodResponse(response):
    while response != "y" and response != "n":
        print("Unkown response")
        response = input("Confirmation? (y/n)")
    return response


def again():
    again = input("Another task? (y/n)")
    while again != "y" and again != "n":
        print("Unknown response")
        again = input("Another task? (y/n)")
    if again == "y":
        # mainMenu()
        print("main menu not faite encore hihi")
        exit()
    else:
        exit()


def waiting():
    """I have to ameliorate this by checking in wich os system this program is running"""
    os.system("cls")
    print("Please open laminate data.xlsx and enter your configurations")
    waiting = input("done? (y/n)")
    while waiting != "y" and waiting != "n":
        print("Unknown response")
        print("Please open laminate data.xlsx and enter your configurations")
        waiting = input("done? (y/n)")
    return waiting


if __name__ == "__main__":
    while True:
        response_fs = forcingGoodResponse(firstScreen())
        while response_fs == "y":
            response_ss = forcingGoodResponse(secondScreen())
            if response_ss == "y":
                thirdScreen()
                exit()
            else:
                wait = waiting()
                if wait == "y":
                    pass
                else:
                    print("Program stopped")
                    exit()
        waiting = waiting()
        if waiting == "y":
            pass
        else:
            print("Programmed stopped")
            exit()
