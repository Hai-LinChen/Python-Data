from math import *
import numpy as np

"""
Batch calculation of molar volume of a phase from its lattice parameters.
    ^ input from and output to files
    ^ using 'dictionary' instead of 'if, elif, else'
    ^ using casefold()
"""

cryslsystems = {'c': 'Cubic', 't': 'tetragonal', 'o': 'orthorhombic', 'r': 'rhombohedral', 'h' : 'hexagonal',
                'm': 'monoclinic', 'tri': 'triclinic'}

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
    elif id == 'c':
        a = float(line[2])
        b, c = a, a
    elif id == 't':
        a = float(line[2])
        c = float(line[3])
        b = a
    elif id == 'o':
        a = float(line[2])
        b = float(line[3])
        c = float(line[4])
    elif id == 'r':
        a = float(line[2])
        b, c = a, a
        alpha = float(line[3])
        beta = alpha
        gamma = alpha
    elif id == 'h':
        a = float(line[2])
        b = a
        c = float(line[3])
        beta = 120
    elif id == 'm':
        a = float(line[2])
        b = float(line[3])
        c = float(line[4])
        beta = float(line[5])
    else:
        print("The data at line {} are wrong!".format(line_number))
    lattice_params = (a, b, c, alpha, beta, gamma)
    n_a_cell = int(line[0])
    n_a_model = int(line[1])
    return n_a_cell, n_a_model, lattice_params

def identification(id): # to identify the crystal structure
    try:
        return cryslsystems[id]
    except Exception:
        return cryslsystems['tri']

def separation(w):  # to sepate the dataset from phase name for each line
    length = len(w)
    if "#" in w:
        idx = w.index("#")
        data_line = w[1:idx]
        if length == idx + 1:
            phasename = 'empty'  # in case "#" is there but nothing exists after it
        else:
            phasename = ''
            for j in range(idx + 1, length): # in case there are 2 or more words
                phasename += w[j]
    else:
        phasename = None
        for i in range(4, length): # in case no space between "#" and the phase name
            if "#" in w[i]:
                data_line = w[1:i]
                phasename = w[i].replace("#", '') # get the phase name (?) from the string "#?"
                try:                              # in case 2 or more words exist after "#"
                    for j in range(i + 1, length):
                        phasename += w[j].replace("#", '')
                except Exception:
                    pass
                break
        if not phasename: # if no "#" has been found
            phasename= 'unspecified'
            data_line = w[1:]
    #print(data_line) # to confirm if this function works.
    return phasename, data_line

def process_data(file):
    f1 = open(file, 'r').readlines()
    f2 = open('volume_ver7.dat', 'w')
    phase = []
    for i in range(len(f1)):
        w=f1[i].split()
        id = w[0].casefold()
        phasename, line = separation(w)
        phase.append(phasename)
        results = analyze_data(id, line, i+1)
        n_a_cell = results[0]
        n_a_model = results[1]
        lattice_params = results[2]
        identity = identification(id)
        molar_volume = round(volume_calculator(lattice_params, n_a_model, n_a_cell), 10)
        f2.write(identity + ' ' + str(n_a_cell) + ' ' + str(n_a_model) + ' ')
        for j in range(len(lattice_params)):
            f2.write(str(lattice_params[j]) + ' ')
        f2.write(str(molar_volume)+'\n')
        print("{:>12}: ".format(phase[i]), end=' ')
        print("{:>12}: ".format(identity), end=' ')
        for i in range(3):
            print("\t{:.4f}".format(lattice_params[i]), end = ' ')
        for i in range(3, 6):
            print("{:>4}".format(round(lattice_params[i])), end = ' ')
        print("\t{:.4E}".format(molar_volume), '\n')

def main():
    process_data('lattice_data_for-ver6.txt')

main()