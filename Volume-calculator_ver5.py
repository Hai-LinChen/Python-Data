from math import *
import numpy as np

"""
Batch calculation of molar volume of a phase from its lattice parameters.
input from and output to files
"""

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

def analyze_data(id, line, line_number):
    alpha, beta, gamma = 90, 90, 90  # enter the default value for the angles
    if(len(line)) == 8:
        a, b, c = float(line[2]), float(line[3]), float(line[4])
        alpha, beta, gamma = float(line[5]), float(line[6]), float(line[7])
    elif id == 'c' or id == 'C':
        a = float(line[2])
        b, c = a, a
    elif id == 't' or id == 'T':
        a = float(line[2])
        c = float(line[3])
        b = a
    elif id == 'o' or id == 'O':
        a = float(line[2])
        b = float(line[3])
        c = float(line[4])
    elif id == 'r' or id == 'R':
        a = float(line[2])
        b, c = a, a
        alpha = float(line[3])
        beta = alpha
        gamma = alpha
    elif id == 'h' or id == 'H':
        a = float(line[2])
        b = a
        c = float(line[3])
        beta = 120
    elif id == 'm' or id == 'M':
        a = float(line[2])
        b = float(line[3])
        c = float(line[4])
        beta = float(line[5])
    else:
        print("The data at line {} are wrong!".format(line_number))
    lattice_params = (a, b, c, alpha, beta, gamma)
    n_a_cell = int(line[0])
    n_a_model = int(line[1])
    print(id, lattice_params, end = ' ')
    return n_a_cell, n_a_model, lattice_params

def load_data(file):
    f1 = open(file, 'r').readlines()
    f2 = open('volume.dat', 'w')
    for i in range(len(f1)):
        w=f1[i].split()
        id = w[0]
        line = w[1:9]
        results = analyze_data(id, line, i+1)
        n_a_cell = results[0]
        n_a_model = results[1]
        lattice_params = results[2]
        molar_volume = round(volume_calculator(lattice_params, n_a_model, n_a_cell), 10)
        calc_volume.append(molar_volume)
        f2.write(id + ' ' + str(n_a_cell) + ' ' + str(n_a_model) + ' ')
        for j in range(len(lattice_params)):
            f2.write(str(lattice_params[j]) + ' ')
        f2.write(str(molar_volume)+'\n')
        print(molar_volume, '\n')

def main():
    global calc_volume
    calc_volume = []
    load_data('lattice_data.txt')

main()