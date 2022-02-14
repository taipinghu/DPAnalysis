#!/usr/bin/env python3


# Constants 
Bohr, Hartree = 0.52917721092, 27.2114
Pi = 3.141592653589793

# Unit convert
bohr2ang = 0.529177249
au2debye = 2.5417462
au2kcal=627.51
au2KJ=2625.5
au2eV=27.2113838
amu2kg=1.66053878e-27
b2m=0.529177249e-10
au2J=4.35974434e-18

Multiwfn = '/home/taipinghu/software/Multiwfn_3.8_dev_src_Linux/Multiwfn'

# Atomic colvence radius
"""
  This header file contains the covalent radii of the atoms in Angstroms from:
  "Covalent radii revisited"  Dalton Trans., 2008, 2832, by
  B. Cordero  V. Gomez, A. E. Platero-Prats, M. Reves, J. Echeverria,
  E. Cremades  F. Barragan and S. Alvarez.  The largest values have
  been chosen for C (sp3) and for Mn  Fe, and Co (high-spin).
- RAK  May 2008"""
#define LAST_COV_RADII_INDEX 96
R_Covalence = (
    2.0,  # ghost?
    0.31,  #H
    0.28,  #He
    1.28,  #Li
    0.96,  #Be
    0.84,  #B
    0.76,  #C
    0.71,  #N
    0.66,  #O
    0.57,  #F
    0.58,  #Ne
    1.66,  #Na
    1.41,  #Mg
    1.21,  #Al
    1.11,  #Si
    1.07,  #P
    1.05,  #S
    1.02,  #Cl
    1.06,  #Ar
    2.03,  #K
    1.76,  #Ca
    1.70,  #Sc
    1.60,  #Ti
    1.53,  #V
    1.39,  #Cr
    1.61,  #Mn
    1.52,  #Fe
    1.50,  #Co
    1.24,  #Ni
    1.32,  #Cu
    1.22,  #Zn
    1.22,  #Ga
    1.20,  #Ge
    1.19,  #As
    1.20,  #Se
    1.20,  #Br
    1.16,  #Kr
    2.20,  #Rb
    1.95,  #Sr
    1.90,  #Y
    1.75,  #Zr
    1.64,  #Nb
    1.54,  #Mo
    1.47,  #Tc
    1.46,  #Ru
    1.42,  #Rh
    1.39,  #Pd
    1.45,  #Ag
    1.44,  #Cd
    1.42,  #In
    1.39,  #Sn
    1.39,  #Sb
    1.38,  #Te
    1.39,  #I
    1.40,  #Xe
    2.44,  #Cs
    2.15,  #Ba
    2.07,  #La
    2.04,  #Ce
    2.03,  #Pr
    2.01,  #Nd
    1.99,  #Pm
    1.98,  #Sm
    1.98,  #Eu
    1.96,  #Gd
    1.94,  #Tb
    1.92,  #Dy
    1.92,  #Ho
    1.89,  #Er
    1.90,  #Tm
    1.87,  #Yb
    1.87,  #Lu
    1.75,  #Hf
    1.70,  #Ta
    1.62,  #W
    1.51,  #Re
    1.44,  #Os
    1.41,  #Ir
    1.36,  #Pt
    1.36,  #Au
    1.32,  #Hg
    1.45,  #Tl
    1.46,  #Pb
    1.48,  #Bi
    1.40,  #Po
    1.50,  #At
    1.50,  #Rn
    2.60,  #Fr
    2.21,  #Ra
    2.15,  #Ac
    2.06,  #Th
    2.00,  #Pa
    1.96,  #U
    1.90,  #Np
    1.87,  #Pu
    1.80,  #Am
    1.69  #Cm
)

#extract from Amesp
atomic_data = [
    [0, 'X', 'X', 0.0], 
    [1, 'H', 'Hydrogen', 1.00782504], 
    [2, 'He', 'Helium', 4.00260324], 
    [3, 'Li', 'Lithium', 7.016003], 
    [4, 'Be', 'Beryllium', 9.0121822], 
    [5, 'B', 'Boron', 11.0093054], 
    [6, 'C', 'Carbon', 12.0], 
    [7, 'N', 'Nitrogen', 14.003074], 
    [8, 'O', 'Oxygen', 15.9949146], 
    [9, 'F', 'Fluorine', 18.9984032], 
    [10, 'Ne', 'Neon', 19.9924356], 
    [11, 'Na', 'Sodium', 22.9897677], 
    [12, 'Mg', 'Magnesium', 23.9850423], 
    [13, 'Al', 'Aluminium', 26.9815386], 
    [14, 'Si', 'Silicon', 27.9769271], 
    [15, 'P', 'Phosphorus', 30.973762], 
    [16, 'S', 'Sulfur', 31.9720707], 
    [17, 'Cl', 'Chlorine', 34.9688527], 
    [18, 'Ar', 'Argon', 39.9623837], 
    [19, 'K', 'Potassium', 38.9637074], 
    [20, 'Ca', 'Calcium', 39.9625906], 
    [21, 'Sc', 'Scandium', 44.95591], 
    [22, 'Ti', 'Titanium', 47.9479473], 
    [23, 'V', 'Vanadium', 50.9439617], 
    [24, 'Cr', 'Chromium', 51.9405098], 
    [25, 'Mn', 'Manganese', 54.9380471], 
    [26, 'Fe', 'Iron', 55.9349393], 
    [27, 'Co', 'Cobalt', 58.9331976], 
    [28, 'Ni', 'Nickel', 57.9353462], 
    [29, 'Cu', 'Copper', 62.9295989], 
    [30, 'Zn', 'Zinc', 63.9291448], 
    [31, 'Ga', 'Gallium', 68.92558], 
    [32, 'Ge', 'Germanium', 73.9211774], 
    [33, 'As', 'Arsenic', 74.9215942], 
    [34, 'Se', 'Selenium', 79.9165196], 
    [35, 'Br', 'Bromine', 78.9183361], 
    [36, 'Kr', 'Krypton', 83.911507], 
    [37, 'Rb', 'Rubidium', 84.911794], 
    [38, 'Sr', 'Strontium', 87.9056188], 
    [39, 'Y', 'Yttrium', 88.905849], 
    [40, 'Zr', 'Zirconium', 89.9047026], 
    [41, 'Nb', 'Niobium', 92.9063772], 
    [42, 'Mo', 'Molybdenum', 97.9054073], 
    [43, 'Tc', 'Technetium', 97.907215],
    [44, 'Ru', 'Ruthenium', 101.904349], 
    [45, 'Rh', 'Rhodium', 102.9055], 
    [46, 'Pd', 'Palladium', 105.903478], 
    [47, 'Ag', 'Silver', 106.905092], 
    [48, 'Cd', 'Cadmium', 113.903357], 
    [49, 'In', 'Indium', 114.903882], 
    [50, 'Sn', 'Tin', 119.902199], 
    [51, 'Sb', 'Antimony', 120.903821], 
    [52, 'Te', 'Tellurium', 129.906229], 
    [53, 'I', 'Iodine', 126.904473], 
    [54, 'Xe', 'Xenon', 131.904144], 
    [55, 'Cs', 'Caesium', 132.905429], 
    [56, 'Ba', 'Barium', 137.905232], 
    [57, 'La', 'Lanthanum', 138.906347], 
    [58, 'Ce', 'Cerium', 139.905433], 
    [59, 'Pr', 'Praseodymium', 140.907647], 
    [60, 'Nd', 'Neodymium', 141.907719], 
    [61, 'Pm', 'Promethium', 144.912743], 
    [62, 'Sm', 'Samarium', 151.919728], 
    [63, 'Eu', 'Europium', 152.921225], 
    [64, 'Gd', 'Gadolinium', 157.924019], 
    [65, 'Tb', 'Terbium', 158.925342], 
    [66, 'Dy', 'Dysprosium', 163.929171], 
    [67, 'Ho', 'Holmium', 164.930319], 
    [68, 'Er', 'Erbium', 165.93029], 
    [69, 'Tm', 'Thulium', 168.934212], 
    [70, 'Yb', 'Ytterbium', 173.938859], 
    [71, 'Lu', 'Lutetium', 174.94077], 
    [72, 'Hf', 'Hafnium', 179.946546], 
    [73, 'Ta', 'Tantalum', 180.947462], 
    [74, 'W', 'Tungsten', 183.950928], 
    [75, 'Re', 'Rhenium', 186.955744], 
    [76, 'Os', 'Osmium', 191.961467], 
    [77, 'Ir', 'Iridium', 192.962917], 
    [78, 'Pt', 'Platinum', 194.964766], 
    [79, 'Au', 'Gold', 196.966543], 
    [80, 'Hg', 'Mercury', 201.970617], 
    [81, 'Tl', 'Thallium', 204.974401], 
    [82, 'Pb', 'Lead', 207.976627], 
    [83, 'Bi', 'Bismuth', 208.980374], 
    [84, 'Po', 'Polonium', 208.982404], 
    [85, 'At', 'Astatine', 209.987126], 
    [86, 'Rn', 'Radon', 222.017571], 
    [87, 'Fr', 'Francium', 223.019733], 
    [88, 'Ra', 'Radium', 226.025403], 
    [89, 'Ac', 'Actinium', 227.02775], 
    [90, 'Th', 'Thorium', 232.038051], 
    [91, 'Pa', 'Protactinium', 231.03588], 
    [92, 'U', 'Uranium', 238.050785], 
    [93, 'Np', 'Neptunium', 237.048168], 
    [94, 'Pu', 'Plutonium', 244.064199], 
    [95, 'Am', 'Americium', 243.061373],    
    [96, 'Cm', 'Curium', 247.070347], 
    [97, 'Bk', 'Berkelium', 247.0703], 
    [98, 'Cf', 'Californium', 251.07958], 
    [99, 'Es', 'Einsteinium', 252.082944], 
    [100, 'Fm', 'Fermium', 257.095099], 
    [101, 'Md', 'Mendelevium', 258.09857], 
    [102, 'No', 'Nobelium', 259.100931], 
    [103, 'Lr', 'Lawrencium', 260.10532], 
    [104, 'Rf', 'Rutherfordium', None], 
    [105, 'Db', 'Dubnium', None], 
    [106, 'Sg', 'Seaborgium', None], 
    [107, 'Bh', 'Bohrium', None], 
    [108, 'Hs', 'Hassium', None], 
    [109, 'Mt', 'Meitnerium', None], 
    [110, 'Ds', 'Darmstadtium', None], 
    [111, 'Rg', 'Roentgenium', None], 
    [112, 'Cn', 'Copernicium', None], 
    [113, 'Uut', 'Ununtrium', None], 
    [114, 'Uuq', 'Ununquadium', None], 
    [115, 'Uup', 'Ununpentium', None], 
    [116, 'Uuh', 'Ununhexium', None], 
    [117, 'Uus', 'Ununseptium', None], 
    [118, 'Uuo', 'Ununoctium', None]
]
name2index = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Uut": 113,
    "Uuq": 114,
    "Uup": 115,
    "Uuh": 116,
    "Uus": 117,
    "Uuo": 118,
}
# END: from phonopy

index2name = {value:key for key, value in name2index.items()}
###########################################################
