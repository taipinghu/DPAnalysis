#!/usr/bin/env python3

import subprocess as sp
import os
import numpy as np
import time

from ase import Atoms

from utils import string_2_array, file_to_list, box_to_cell, str_to_type_map

###########################################################
def fileio_functions():
    while True:
        print("\n")
        print("========================================================")
        print("                    Main Menu 10                        ")
        print("                  File I/O Options                      ")
        print(" (  1)  Traj combination of LAMMPS dump file            ")
        print(" (  2)  Convert a multi-frames .dump file to .pdb file  ")
        print("                                                        ")
        print("               Other format to LAMMPS                   ")
        print(" ( 10)  Convert POSCAR to LAMMPS lmp file              ")
        print(" ( 99)  Convert POSCAR to standard .car format file     ")
        print(" (100)  Convert multiple POSCAR to traj.xyz file        ")
        print("\n Tips: Input -10 to return main menu                  ")
        print("========================================================")
        jsel = input("\n Please input the menu index:\n")
    
        if jsel == str(-10):
            break
        elif jsel == str(1):
            str_input = input("\n Please input the begin frame, end frame and step number, e.g., 0, 10, 1\n")
            begin, end, step = string_2_array(str_input)
            print(begin, end, step)
            for iframe in range(begin, end, step):
                assert(os.path.exists('%d.lammpstrj'))
                sp.call('cat '+ '%d.lammpstrj'%iframe + ' >> test.dump', shell=True)
                print(' %d.lammpstrj has been exported...'%iframe)
        elif jsel == str(2):
            pass
        elif jsel == str(10):
            str_in = input("\n Please input the type map, e.g., Li, Si. ")
            convert_poscar_to_lmp(type_map=str_to_type_map(str_in))
        else:
            print(" !!! Unspported menu index, please input again. !!!")
###########################################################


###########################################################
def POSCAR_to_car(file_name="POSCAR", out_file="test.car"):
    if not os.path.exists(file_name):
        raise ValueError(" The file_name is not exists, please check again.")
    sp.call("ase -T convert " + file_name + " " + out_file, shell=True)
###########################################################

###########################################################
def convert_poscar_to_lmp(type_map):
    import dpdata

    system = dpdata.System("POSCAR", fmt="vasp/poscar", type_map=type_map)
    system.to_lammps_lmp("POSCAR.lmp")
###########################################################
