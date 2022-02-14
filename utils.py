#!/usr/bin/env python3

import glob
import numpy as np
import shutil
import os
from collections import Counter
import logging
import sys
from functools import lru_cache as cache
from copy import copy
from itertools import count

import define

SHORT_CMD="DPAnalysis"
DPAlog = logging.getLogger(__name__)
DPAlog.setLevel(logging.INFO)
DPAlogf = logging.FileHandler(os.getcwd()+os.sep+SHORT_CMD+'.log', delay=True)
DPAlogf_formatter=logging.Formatter('%(asctime)s - %(levelname)s : %(message)s')
DPAlogf.setFormatter(DPAlogf_formatter)
DPAlog.addHandler(DPAlogf)

###########################################################
class Box(object):
    def __init__(self, box, moles):
        self.__box = box
        self.__moles = moles       

    def get_nmoles(self):
        return len(self.__moles)

    def genCluster(self, centers, radius):
        @cache(maxsize=None)
        def genCellVectors(n):
            cells = set()
            for i in range(0, n + 1):
                for j in range(0, n + 1 - i):
                    k = n - i - j
                    cells.add(( i,  j,  k)); cells.add((-i, -j, -k))
                    cells.add((-i,  j,  k)); cells.add(( i, -j, -k))
                    cells.add(( i, -j,  k)); cells.add((-i,  j, -k))
                    cells.add(( i,  j, -k)); cells.add((-i, -j,  k))
            return np.einsum('ij,jk->ik', np.array([cell for cell in cells]), self.__box)

        nmole = self.get_nmoles()
        moles = [self.__moles[center].move(np.array([0, 0, 0])) for center in centers]
        newcenters = list(range(len(centers)))
        masscenters = np.array([mole.getMassCenter() for mole in self.__moles], dtype=np.double)

        if self.__box is None:
            pass
        else:
            for mole_index in range(nmole):
                notfound = 0
                curMole = self.__moles[mole_index]
                for i in count():
                    notfound += 1
                    rn = []
                    aimCellVector = genCellVectors(i)
                    for center in centers:
                        r = aimCellVector + masscenters[mole_index] - masscenters[center]
                        rn.append(np.sqrt(np.einsum('ij,ij->i', r, r)))
                    rn = np.min(rn, axis=0)
                    for index in np.where((rn <= radius) & (rn > 0.01))[0]:
                        moles.append(curMole.move(aimCellVector[index]))
                        notfound = 0
                    if notfound > 3: break
        return len(moles)


class Molecule(object):
    def __init__(self, Symbols, Coordinates, Charges, ResName=None, Energy=0.0, AtomicMMTypes=None, Connectivities=None):
        self.__AtomicSymbols = tuple(Symbols)
        self.__AtomicNumbers = tuple([define.name2index[symbol] for symbol in self.__AtomicSymbols])
        self.__AtomicMasses = tuple([define.atomic_data[num][3] for num in self.__AtomicNumbers])
        self.__AtomicCharges = tuple(Charges)
        self.__AtomicCoordinates = np.array(Coordinates, dtype=np.double)

        if ResName is not None: self.__ResName = ResName
        self.__Energy = Energy
        if AtomicMMTypes is not None: 
            self.__AtomicMMTypes = tuple(AtomicMMTypes)
        else:
            self.__AtomicMMTypes = []
        if Connectivities is not None: 
            self.__MoleConnectivities = tuple(Connectivities)   
        else:
            self.__MoleConnectivities = []

        self.__MoleCharge = 0
        self.__MoleSpin = 0.0

        weights = np.array(self.__AtomicMasses)
        self.__MoleMass = weights.sum()
        self.__MassCenter = np.einsum('ij,i->j', self.__AtomicCoordinates, weights) / self.__MoleMass

    def getAtomicMasses(self):
        return self.__AtomicMasses

    def setAtomicCharges(self, list):
        self.__AtomicCharges = list
        return self.__AtomicCharges

    def setResName(self, string):
        self.__ResName = string
        return self.__ResName

    def getNAtoms(self):
        return len(self.__AtomicSymbols)
    
    def getResName(self):
        return self.__ResName

    def getMoleConnectivities(self):
        return self.__MoleConnectivities

    def getMoleCharge(self):
        return self.__MoleCharge

    def getMoleSpin(self):
        return self.__MoleSpin

    def setMoleCharge(self, charge):
        self.__MoleCharge = charge

    def setMoleSpin(self, spin):
        self.__MoleSpin = spin

    def addMoleCharge(self, charge):
        self.__MoleCharge += charge
        chgperatom = charge / self.natoms()
        self.__AtomicCharges = tuple([c + chgperatom for c in self.__AtomicCharges])
        return self.getMoleCharge()

    def getAtomicCoordinates(self, move=np.array([0.0, 0.0, 0.0], dtype=np.double), factor=1):
        return (self.__AtomicCoordinates + move) * factor

    def getAtomicSymbols(self):
        return self.__AtomicSymbols

    def getAtomicNumbers(self):
        return self.__AtomicNumbers

    def getAtomicCharges(self):
        return self.__AtomicCharges

    def getAtomicMMTypes(self):
        return self.__AtomicMMTypes

    def getMoleMass(self):
        return self.__MoleMass

    def getMoleEnergy(self):
        return self.__Energy

    def getMassCenter(self):
        return self.__MassCenter.copy()

    def renewAtomicCharges(self, newCharges):
        self.__AtomicCharges = tuple(newCharges)

    def renewAtomicCoordinates(self, coordinates):
        self.__AtomicCoordinates = np.array(coordinates, dtype=np.double)

    def move(self, vector=np.array([0.0, 0.0, 0.0], dtype=np.double)):
        ret = copy(self)
        ret.__AtomicCoordinates = ret.__AtomicCoordinates + vector
        ret.__MassCenter = ret.__MassCenter + vector
        return ret

###########################################################


###########################################################
def get_time_fp(filename):
    """
    This code is written by dfz

    Adopted at 20211217
    """   
    with open(filename, "r") as fip:
        lines = fip.readlines()
    for ii in lines[::-1]:
        if "Total CPU time used (sec):" in ii:
            return float(ii.split()[5])/3600
###########################################################


###########################################################
def prange(start, end, step):
    if start < end:
        for i in range(start, end, step):
            yield i, min(i+step, end)
###########################################################


###########################################################
def process_bar(I, Max, Process=' Reading File:'):
    percent = I * 100.0 / Max
    e = int(I * 20 / Max)
    y = 20 - e
    process_bar = '{:s}'.format(Process) + '[' + '>' * e + '-' * y + ']' + '{:.2f}'.format(percent) + '%' +  '\r'
    sys.stdout.write(process_bar)
    sys.stdout.flush()
###########################################################


###########################################################
def symbols_to_typemap(symbols: list) -> dict:
    type_map = {}
    for sym in symbols:
        type_map[sym] = type_map.get(sym, 0) + 1
    return type_map
###########################################################


###########################################################
def typemap_list_to_symbols(atom_numbs: list, atom_names: list) -> list:
    atomic_symbols = []
    idx = 0
    for numb in atom_numbs:
        atomic_symbols.extend((atom_names[idx], )*numb)
        idx += 1
    return atomic_symbols
###########################################################


###########################################################
def typemap_to_symbols(Dict: dict) -> list:
    if not isinstance(Dict, dict):
        raise TypeError(' The argument must be the type of dict.')
    symbols = []
    _atom_names = list(Dict.keys())
    _atom_numbs = list(Dict.values())
    idx = 0
    for numb in _atom_numbs:
        symbols.extend((_atom_names[idx], )*numb)
        idx += 1
    return symbols
###########################################################


###########################################################
def box_to_cell(box):
    """
    Convert Box to Cell a, b, c, alpha, beta, gamma
    """
    if not isinstance(box, np.ndarray):
        raise TypeError('The argument Box must be the type of np.array')

    a = np.sqrt(box[0,0]**2+box[0,1]**2+box[0,2]**2)
    b = np.sqrt(box[1,0]**2+box[1,1]**2+box[1,2]**2)
    c = np.sqrt(box[2,0]**2+box[2,1]**2+box[2,2]**2)

    alpha = np.arccos((box[0,0]*box[1,0]+box[0,1]*box[1,1]+box[0,2]*box[1,2])/a/b) * 180 / np.pi
    beta  = np.arccos((box[0,0]*box[2,0]+box[0,1]*box[2,1]+box[0,2]*box[2,2])/a/c) * 180 / np.pi
    gamma = np.arccos((box[1,0]*box[2,0]+box[1,1]*box[2,1]+box[1,2]*box[2,2])/b/c) * 180 / np.pi
    
    return a, b, c, alpha, beta, gamma
###########################################################


###########################################################
def select_dirs(root_path: str) -> list:
    """
    Return all the sub-dirs in root_path

    Written at 20211225
    """
    select_dirs = []
    for tmp in glob.glob(os.path.join(root_path, "*")):
        if os.path.isdir(tmp):
            select_dirs.append(tmp)
    return select_dirs
###########################################################


###########################################################
def select_all_dirs(root_path: str, dir_name: str) -> list:
    select_dirs = []
    for r, d, fs in os.walk(root_path):
        for tmp_dir in d:
            if os.path.islink(tmp_dir):
                continue
            if tmp_dir == dir_name:
                select_dirs.append(os.path.join(r, tmp_dir))
    return select_dirs
###########################################################


###########################################################
def select_logs(root_path, fname='OUTCAR'):
    """
    Find all VASP OUTCAR file in root_path and its all sub-dirs
    This code is written by dfz

    Adopted at 20211217
    """
    logs = []
    for r, d, fs in os.walk(root_path):
        for f in fs:
            if os.path.islink(f):
                continue
            if f == fname:
                logs.append(os.path.join(r, f))
    return sorted(logs)
###########################################################


###########################################################
def get_file_line_number(file_name):
    if not os.path.getsize(file_name):
        return 0
    counter = -1
    for counter, line in enumerate(open(file_name, 'rU')):
        pass
        counter += 1
    return counter
###########################################################


###########################################################
def figure_output_options(plt):
    while True:
        ksel = figure_utils()
        if ksel == 0:
            plt.show()
        elif ksel == 1:
            plt.savefig(exported_figure_file())
            print(" Export finished!")
            break
        elif ksel == -10:
            break
    return 
##########################################################


###########################################################
def figure_utils():
    print(" (0)  Display figure in screen")
    print(" (1)  Save figure in current dir")
    print("\n Tips: Input -10 to return main menu")
    ksel = int(input().strip())
    return ksel
###########################################################


###########################################################
def exported_figure_file():
    file_export = input(" Please input the filename to export:\n").strip()
    return file_export    
###########################################################


###########################################################
def create_path (path) :
    """
    Create a path, adapted from dpgen code
    """
    path += '/'
    if os.path.isdir(path) : 
        dirname = os.path.dirname(path)        
        counter = 0
        while True :
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname) : 
                shutil.move (dirname, bk_dirname) 
                break
                counter += 1
    os.makedirs (path)

###########################################################


###########################################################
def string_2_array(inpstr):
    '''
    Parse input integer string (e.g. 3,4,5,7-29,99) to array
    
    e.g. inputstring: 3,4,5,7-9,99
         return: 3,4,5,7,8,9,99
    
    Written at 20190921
    '''
    array = []
    tmp1 = str(inpstr).split(',')
    for i in range(len(tmp1)):
        if '-' in tmp1[i]:
            tmp2 = tmp1[i].split('-')
            for j in range(int(tmp2[0]), int(tmp2[1])+1):
                array.append(j)
        else:
            array.append(int(tmp1[i]))
    length = len(array)
    return array
###########################################################


###########################################################
def KDE_Func(xdat, ydat, xstep=0.02):
    minval, maxval = np.min(xdat), np.max(xdat)
    num = round((maxval - minval) / xstep) + 1
    alpha = 1 / (2 * 0.5**2)

    factor = 1.0 / len(xdat) * np.sqrt(alpha / np.pi)
    x = np.linspace(minval, maxval, int(num))
    result = np.zeros((len(x), 4))
    result[:, 0] = x[:]
    t = np.zeros(xdat.shape)
    for idx, xval in enumerate(x):
        np.exp((xdat - xval)**2 * (-alpha), out=t)
        result[idx, 1] = factor * np.sum(t)
        result[idx, 2] = factor * np.sum(t * ydat) / result[idx, 1]
        result[idx, 3] = np.sqrt(factor * np.sum((ydat - result[idx, 2])**2 * t) / result[idx, 1])
    return result[:, (0, 1, 2, 3)]
###########################################################


###########################################################
def local_label(file_name, label, mode=1):
    '''
    Locate the line where the label first appears.

    Arguments:
        file_name:  str     file name
        label:      str     label to be found
        mode:       int     =1: return the line, =2: return the nline

    Written at 20190901
    '''
    with open(file_name, 'r') as f:
        for nline, line in enumerate(f):
            if label in line:
                if mode == 1:
                    return line
                if mode == 2:
                    return nline
###########################################################


###########################################################
def file_to_list(file_name):
    with open(file_name, 'r') as f:
        lines = [line.rstrip() for line in f.readlines()]

    return lines
###########################################################


###########################################################

###########################################################


###########################################################
def box_to_cell(box):
    """
    Convert Box to Cell a, b, c, alpha, beta, gamma
    """
    if not isinstance(box, np.ndarray):
        raise TypeError('The argument Box must be the type of np.array')

    a = np.sqrt(box[0,0]**2+box[0,1]**2+box[0,2]**2)
    b = np.sqrt(box[1,0]**2+box[1,1]**2+box[1,2]**2)
    c = np.sqrt(box[2,0]**2+box[2,1]**2+box[2,2]**2)

    alpha = np.arccos((box[0,0]*box[1,0]+box[0,1]*box[1,1]+box[0,2]*box[1,2])/a/b) * 180 / np.pi
    beta  = np.arccos((box[0,0]*box[2,0]+box[0,1]*box[2,1]+box[0,2]*box[2,2])/a/c) * 180 / np.pi
    gamma = np.arccos((box[1,0]*box[2,0]+box[1,1]*box[2,1]+box[1,2]*box[2,2])/b/c) * 180 / np.pi
    
    return a, b, c, alpha, beta, gamma
###########################################################


###########################################################
def cell_to_box(a, b, c, alpha, beta, gamma):
    alpha = alpha / 180 * np.pi
    beta  = beta  / 180 * np.pi
    gamma = gamma / 180 * np.pi

    box = np.zeros((3,3), dtype=np.double) 
    box[0, 0] = a
    box[0, 1] = 0
    box[0, 2] = 0
    box[1, 0] = b * np.cos(gamma)
    box[1, 1] = b * np.sin(gamma)
    box[1, 2] = 0
    box[2, 0] = c * np.cos(beta)
    box[2, 1] = c * (np.cos(alpha)-np.cos(beta)*np.cos(gamma)) / np.sin(gamma)
    box[2, 2] = c * np.sqrt(1 - np.cos(beta)**2 - ((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma))**2)
    return box
###########################################################


###########################################################
def box_to_recvec(box):
    a1 = box[0]
    a2 = box[1]
    a3 = box[2]
    b1 = np.cross( a2, a3 )
    b2 = np.cross( a3, a1 )
    b3 = np.cross( a1, a2 )
    volume = np.dot( a1, np.cross( a2, a3 ) )
    rec_vec = [ b1, b2, b3 ]
    # it follows the definition for b_j: a_i * b_j = 2pi * delta(i,j)
    for i in range(0,3):
        for j in range(0,3):
            rec_vec[i][j] = rec_vec[i][j] * 2 * np.pi / volume
    return rec_vec
###########################################################


###########################################################
def typemap_to_symbols(type_map):
    """
    Convert a Type map dict to symbol list.
    
    Written at 20210909
    """
    if not isinstance(type_map, dict):
        raise TypeError(' The argument must be the type of dict.')
    symbols = []
    _atom_names = list(type_map.keys())
    _atom_numbs = list(type_map.values())

    idx = 0
    for numb in _atom_numbs:
        symbols.extend((_atom_names[idx], )*numb)
        idx += 1
    return symbols
###########################################################


########################################################### 
def arg_sort(iterable, reverse=False):
    length = len(iterable)
    return sorted(range(length), key=lambda index: round(iterable[index], 8), reverse=reverse)
###########################################################


###########################################################
def dumpbox_to_box(dumpbox):
    bounds = np.zeros([3,2])
    tilt = np.zeros([3])
    for dd in range(3):
        info = dumpbox[dd]
        bounds[dd][0] = info[0]
        bounds[dd][1] = info[1]
        tilt[dd] = info[2]

    xy, xz, yz = tilt[0], tilt[1], tilt[2]
    xlo = bounds[0][0] - min(0.0,xy,xz,xy+xz)
    xhi = bounds[0][1] - max(0.0,xy,xz,xy+xz)
    ylo = bounds[1][0] - min(0.0,yz)
    yhi = bounds[1][1] - max(0.0,yz)
    zlo = bounds[2][0]
    zhi = bounds[2][1]
    info = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]

    orig = np.array([info[0][0], info[1][0], info[2][0]])
    lens = []
    for dd in range(3) :
        lens.append(info[dd][1] - info[dd][0])
    xx = [lens[0], 0, 0]
    yy = [xy, lens[1], 0]
    zz=  [xz, yz, lens[2]]
    return np.array([xx, yy, zz])
###########################################################
    
###########################################################
##             Those programs are used in DP             ##
###########################################################
def make_iter_name(iter_index):
    return "iter." + "%06d"%(iter_index)
###########################################################


###########################################################
##                Other Utils                           ###
###########################################################
def box_to_volume(box):
    """
    Calculate the volume of the box, using formula V = np.dot(a, np.cross(b, c))
    """

    if not isinstance(box, np.ndarray):
        raise TypeError('The argument Box must be the type if np.ndarray')

    return np.dot(box[0], np.cross(box[1], box[2]))

def str_to_type_map(str_in: str) -> list:
    """
    Convert a string into type map list. e.g., input "Li, Si", return ["Li", "Si"]

    Args:
        str_in : input string

    Written at 20220119, by Taiping Hu.
    """
    type_map = [tmp.strip() for tmp in str_in.split(",")]
    return type_map
###########################################################
    



