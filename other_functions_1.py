#!/usr/bin/env python3

from utils import get_time_fp, box_to_volume, box_to_cell, select_logs
from fileio import read_file
from define import amu2kg

def other_functions_1():
    while True:
        print("\n")
        print("=======================================================")
        print("                    Other Functions I                  ")
        print("                  Varisou Analysis Tools               ")
        print(" (  1)  Print time cost in sing VASP task              ")
        print("                                                       ")
        print("                   Cell Informations                   ")
        print(" ( 90)  Print system information                       ")
        print(" ( 91)  Print box info in current all sub-dirs         ")
        print("\n Tips: Input -10 to return main manu                 ")
        print("=======================================================")
        jsel = input("\n Please input the menu index: ")
    
        if jsel == str(-10):
            break
        elif jsel == str(1):
            cost_time = get_time_fp("OUTCAR")
            print(" Time cost: %.2f h"%cost_time)
        elif jsel == str(90):
            file_in = input("\n Please input the file name: ")
            format_in = input("\n Please input the file format: ")
            box, moles = read_file(file_in, format_in)
            print("")
            print(" Box lattice constant: {}".format(box_to_cell(box)))
            print(" Box Volume: %.5f A^3"%box_to_volume(box))
            print(" Volume per atom: %.5f A^3"%(box_to_volume(box)/len(moles)))
            mass = sum([mole.getMoleMass() for mole in moles])
            print(" Total Mass: %.5f Kg"%(mass))
            print(" Density: %.5f g/cm^3"%(mass*amu2kg*1e3/(box_to_volume(box)*1e-24)))
        elif jsel == str(91):
            str_in = input("\n Please input the root dir name: ")
            format_in = input("\n please input the file format (vasp/poscar, vasp/contcar, vasp/outcar):")
            if format_in == "vasp/poscar":
                all_files = select_logs(str_in, "POSCAR")
            elif format_in == "vasp/contcar":
                all_files = select_logs(str_in, "CONTCAR")
            
            if len(all_files) == 0:
                print(" Cannot find the POSCAR file in current dir, please check again.")
                break
            for tmp_file in all_files:
                box, moles = read_file(tmp_file, format_in)
                print(tmp_file, box_to_cell(box))

            

            
