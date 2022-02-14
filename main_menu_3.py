#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt

from utils import select_all_dirs, file_to_list
from plot_new import bold_axis_and_ticks, set_legend

###########################################################
def dpgen_autotest_functions():
    while True:
        print("\n")
        print("========================================================")
        print("                   DPGEN autotest module                ")
        print("                                                        ")
        print("       Post processing functions for dpgen autotest     ")
        print(" (  1)  Perform dpgen autotest (eos and elastic)        ")
        print(" ( 10)  Plot EOS figure                                 ")
        print("\n Tips: Input -10 to return main menu                  ")
        print("========================================================")

        jsel = input("\n Please input the menu index: ")

        if jsel == str(-10):
            break
        elif jsel == str(10):
            str_in = input("\n Please input the root dir name: ")
            plot_eos_figure(str_in)
###########################################################


###########################################################
def plot_eos_figure(root_dir):
    eos_task_total = select_all_dirs(root_dir, "eos_00")
    
    vasp_x, vasp_y = [], []
    dp_x, dp_y = [], []
    for eos_task in eos_task_total:
        result = os.path.join(eos_task, "result.out")
        lines = file_to_list(result)
        for line in lines[2:]:
            tmp = list(map(float, line.split()))
            if "vasp" in eos_task:
                vasp_x.append(tmp[0])
                vasp_y.append(tmp[1])
            elif "deepmd" in eos_task:
                dp_x.append(tmp[0])
                dp_y.append(tmp[1])
    fig = plt.figure()
    axes = fig.add_subplot()
    axes.plot(vasp_x, vasp_y, marker="o", markersize=12)
    axes.plot(dp_x, dp_y, marker="o", markersize=12)
    set_legend(axes, ["DFT", "DP"], loc=1)
    bold_axis_and_ticks(axes)
    plt.savefig("EOS.png")
    #plt.show()
###########################################################