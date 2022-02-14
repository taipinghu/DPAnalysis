#!/usr/bin/env python3

import json
import numpy as np
import os, glob, sys
import warnings
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sb
from matplotlib.colors import LogNorm
import logging

print("####################################################")
print("##                    DPAnalysis                  ##")
print("##           Deep Potential Analysis tools        ##")
print("##                 Version 0.0.1                  ##")
print("##    Developed by Battery Team of DP Technology  ##")
print("##           Email: taipinghu@iccas.ac.cn         ##")
print("####################################################\n")

#first delete the existed DPAnalysis.log file
#if os.path.exists("DPAnalysis.log"):
#    os.remove("DPAnalysis.log")

while True:
    print("\n")
    print("=======================================================")
    print("                   Main Menu                          ")
    print(" (  1)  Various functions used in dpgen iteration     ")
    print(" (  2)  Various functions used in deepmd-kit          ")
    print(" (  3)  DPGEN autotest module                         ")
    print(" (  4)  Deep potential molecular dynamics (DPMD) Analysis")
    print(" (  5)  Various functions used in fp calculations     ")
    print(" (  6)  Various functions processing Box class        ")
    print(" (100)  File I/O Options                              ")
    print(" (300)  Other function I                              ")
    print("\n Tips: Input q or -10 to exit program               ")
    print("=======================================================")
    isel = input("\n Please input the manu index: ")

    if isel == "q" or isel == str(-10):
        sys.exit()
    elif isel == str(1):
        from main_menu_1 import dpgen_iteration_functions
        dpgen_iteration_functions()
    elif isel == str(2):
        from main_menu_2 import dp_test_functions
        dp_test_functions()
    elif isel == str(3):
        from main_menu_3 import dpgen_autotest_functions
        dpgen_autotest_functions()
    elif isel == str(4):
        from main_menu_4 import DPMD_functions
        DPMD_functions()
    elif isel == str(5):
        from main_menu_5 import fp_task_functions
        fp_task_functions()
    elif isel == str(6):
        from main_menu_6 import box_functions
        box_functions()
    elif isel == str(100):
        from main_menu_100 import fileio_functions
        fileio_functions()
    elif isel == str(9):
        from taiping_code import taiping_self_functions
        taiping_self_functions()
    elif isel == str(300):
        from other_functions_1 import other_functions_1
        other_functions_1()
    else:
        print(" !!! Unspported menu index, please input again. !!!")
