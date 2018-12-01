from math import *

"""
Calculation of molar volume of a phase from its lattice parameters.
input from keyboard and output to screen
"""

def enter_params(cry_str):
    alpha, beta, gamma = 90, 90, 90   # enter the default value for the angles
    if cry_str == 'c' or cry_str == 'C':
        a = float(input("enter lattice parameter (in nm): a "))
        b, c = a, a
    elif cry_str == 't' or cry_str == 'T':
        print("enter lattice parameter (in nm): ")
        a = float(input("a: "))
        c = float(input("c: "))
        b = a
    elif cry_str == 'o' or cry_str == 'O':
        print("enter lattice parameter (in nm): ")
        a = float(input("a: "))
        b = float(input("b: "))
        c = float(input("c: "))
    elif cry_str == 'r' or cry_str == 'R':
        a = float(input("enter lattice parameter (in nm): a "))
        b, c = a, a
        alpha = float(input("enter lattice parameter (degree): alpha "))
        beta = alpha
        gamma = alpha
    elif cry_str == 'h' or cry_str == 'H':
        print("enter lattice parameters (in nm): ")
        a = float(input("a: "))
        b = a
        c = float(input("c: "))
        beta = 120
    elif cry_str == 'm' or cry_str == 'M':
        print("enter lattice parameter (in nm): ")
        a = float(input("a: "))
        b = float(input("b: "))
        c = float(input("c: "))
        beta = float(input("enter lattice parameter (degree): beta "))
    else:
        print("enter lattice parameter (in nm): ")
        a = float(input("a: "))
        b = float(input("b: "))
        c = float(input("c: "))
        print("enter lattice parameter (degree): ")
        alpha = float(input("alpha: "))
        beta = float(input("beta: "))
        gamma = float(input("gamma: "))
    lattice_params = (a, b, c, alpha, beta, gamma)
    return lattice_params

def volume_calculator(parameters, na_model, na_cell):
    a, b, c, alpha, beta, gamma = parameters
    alpha = alpha/180*3.14159265
    beta = beta/180*3.14159265
    gamma = gamma/180*3.14159265
    cos_a = cos(alpha)
    cos_b = cos(beta)
    cos_g = cos(gamma)
    cos_a_sq = cos_a * cos_a
    cos_b_sq = cos_b * cos_b
    cos_g_sq = cos_g * cos_g
    volume = (a*b*c) * sqrt(1-cos_a_sq-cos_b_sq-cos_g_sq+2*cos_a*cos_b*cos_g)
    N_avgadro = 6.02214179e+23  # Avgadro number
    molar_volume = (volume * na_model /na_cell) * 1e-27 * N_avgadro
    return molar_volume

def main():
    print("------------------------------------------------------------------")
    print(" c - cubic \n t - tetragonal \n o - orthorhombic \n r - rhombohedral \n h - hexagonal \n m - monoclinic "
          "\n * - triclinic")
    cryst_type = input("Please enter the crystal structure type: ")
    n_a_cell = int(input("Number of atoms in the unit cell: "))
    n_a_model = int(input("Number of atoms in the formula of the model: "))
    lattice_parameters = enter_params(cryst_type)
    molar_volume = volume_calculator(lattice_parameters, n_a_model, n_a_cell)
    print(round(molar_volume, 10))
    print("------------------------------------------------------------------")

main()