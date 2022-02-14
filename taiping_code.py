#!/usr/bin/env python3

#import global module
import os, glob, random, shutil, sys, json
import subprocess as sp
from matplotlib.pyplot import title
import numpy as np
from matplotlib.colors import LogNorm
from numpy.lib.index_tricks import OGridClass
#import my module
from utils import DPAlog, create_path, arg_sort, string_2_array, make_iter_name, select_logs, typemap_to_symbols, symbols_to_typemap
from utils import select_dirs, select_all_dirs

###########################################################
def taiping_self_functions():
    while True:
        print("=============================================================")
        print("       Warning: This module is only used by Taiping Hu      ")
        print("                        Main menu 9                         ")
        print(" (  1) Generate defect structure for specified POSCAR file  ")
        print(" (  2) Generate anti-defect structure for specified POSCAR file")
        print(" ( 97) Structure pertubation                               ")
        print(" ( 98) Generate NEB path and perted structures              ")
        print(" ( 99) dp test using DP calculator                          ")
        print(" (100) Sorted POSCAR according to specified type map        ")
        print(" (999) Plot energy comparision figures based on 99 function ")
        print("\n Tips: Input -10 to return main menu or q to exit program ")
        print("=============================================================")
        jsel = input("\n Please input the menu index:\n")
        if jsel == "q":
            sys.exit()
        if jsel == str(1):
            total_struct_numb = int(input("\n Please input the total numbers of defect structures: "))
            output_dir = input("\n Please input the dir to export files:  ")
            generate_defect(ori_poscar="POSCAR", total_struc_numb=total_struct_numb, output_dir=output_dir)
        elif jsel == str(3):
            pass
        elif jsel == str(93):
            extract_different_AIMD_traj_from_init_data(scale=1.000)
        elif jsel == str(98):
            print("\n This function firstly make NEB path calling vtst scripts,")
            print(" and then perform pertubation for intermediate images.\n")
            json_file = input("\n Please input the parameter file (.json format): ")
            if not os.path.exists(json_file):
                raise FileNotFoundError(" Cannot find %s file, plaese check again."%json_file)
            pertubation_cineb(json_file)

        elif jsel == str(2):
            total_struct_numb = int(input("\n Please input the total numbers of defect structures: "))
            output_dir = input("\n Please input the dir to export files: ")
            generate_anti_defect(ori_poscar="POSCAR", total_struct_numb=total_struct_numb, output_dir=output_dir)
        elif jsel == str(100):
            if not os.path.exists("POSCAR"):
                print("\n Cannot find POSCAR file in current dir, please check again.\n")
                continue
            import dpdata
            sys_dpdata = dpdata.System("POSCAR", fmt="vasp/poscar")
            print(" The old atom name list: {}".format(sys_dpdata.get_atom_names()))
            print(" The old type map dict : {}".format(symbols_to_typemap(sys_dpdata.get_atom_names())))
            import ast
            str_in = ast.literal_eval(input("\n Please   input the new type map, it is a dict, e.g., {\"Li\":10, \"Co\":12, \"O\":24} \n"))
            symbols = typemap_to_symbols(str_in)
            sys_dpdata.sort_atom_names(symbols)
            print("\n The new atom name list: {}".format(sys_dpdata.get_atom_names()))
            print(" The new type map dict : {}".format(symbols_to_typemap(sys_dpdata.get_atom_names())))
            sys_dpdata.to_vasp_poscar("POSCAR_New")

        elif jsel == str(99):
            str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
            for iter_idx in string_2_array(str_in):
                all_tasks = glob.glob(os.path.join(make_iter_name(iter_idx), "02.fp/task.*"))
                if len(all_tasks) == 0:
                    print(" Cannot find task dirs in iter.%06d..\n"%iter_idx)
                    continue
                print(" Calculating in iter.%06d..."%iter_idx)
                dp_test_using_dp_calculator(iter_idx)
        elif jsel == str(999):
            total_json_files = glob.glob(os.path.join("./", "results_iter.*.json"))
            total_json_files.sort()
            total_temps = []
            for json_file in total_json_files:
                jdata = json.load(open(json_file))
                total_temps.extend(jdata.keys())
            total_temps = sorted(set(total_temps))
            
            import matplotlib.pyplot as plt
            fig = plt.figure()
            axes = fig.add_subplot()
            data_DFT = {}
            data_DP  = {}
            for temp in total_temps:
                data_DFT[temp] = []
                data_DP[temp]  = []
                                   
            for json_file in total_json_files:
                jdata = json.load(open(json_file))
                for temp in total_temps:
                    if temp in jdata.keys():
                        data_DFT[temp].extend(jdata[temp]["DFT"])
                        data_DP[temp].extend(jdata[temp]["DP"])
            
            for temp in total_temps:
                axes.scatter(data_DFT[temp], data_DP[temp], alpha=0.6)
            from plot_new import set_legend, bold_axis_and_ticks
            bold_axis_and_ticks(axes)
            set_legend(axes, ["%s K"%temp for temp in total_temps], loc=2)
            axes.set_aspect(1)
            plt.show()
        elif jsel == str(97):
            json_file = input("\n Please input the json file: ")
            structure_pertubation(json_file)
        elif jsel == str(-10):
            break
###########################################################

###########################################################
##   Those scripts were used in CI-NEB and pertubation   ##
###########################################################
def pertubation_cineb(json_file):
    jdata = json.load(open(json_file))

    stage_list = [int(i) for i in jdata['stages']]
    for stage in stage_list:
        if stage == 1:
            print(' Current stage is 1, generate NEB directory...')
            create_neb_dir(jdata)
        elif stage == 2:
            print(' Current stage is 2, randomly perturbate structure...')
            pert_scaled(jdata)
        elif stage == 3:
            print(' Current stage is 3, submit tasks and perform VASP NEB short step (e.g., 20 steps) optimization...')
        elif stage == 4:
            print(' Current stage is 4, collect data from VASP OUTCAR...')
        else:
            raise RuntimeError("unknown stage %d" % stage)

def create_neb_dir(jdata):
    initial_state = jdata['initial_path']
    final_state   = jdata['final_path']
    neb_numb      = jdata['neb_numb']
    out_dir       = jdata['out_dir']
    neb_dir       = jdata['neb_dir']

    cwd = os.getcwd()
    init_path = os.path.join(out_dir, '00.nebmake')
    init_path = os.path.abspath(init_path)
    if not os.path.isdir(init_path):
        create_path(init_path)
    os.chdir(init_path)

    neb_make_cmd = os.path.join(neb_dir, 'nebmake.pl ')
    neb_make_cmd = neb_make_cmd + initial_state + ' ' + final_state + ' ' + str(neb_numb)
    neb_make_cmd = 'perl ' + neb_make_cmd
    sp.check_call(neb_make_cmd, shell=True)
    os.chdir(cwd)

def pert_scaled(jdata):
    out_dir = jdata['out_dir']
    pert_box = jdata['pert_box']
    pert_atom = jdata['pert_atom']
    pert_numb = jdata['pert_numb']
    neb_numb = jdata['neb_numb']

    incar = jdata['incar']
    potcar = jdata['potcar']
    #kpoints = jdata['kpoints']

    init_path = os.path.join(out_dir, '00.nebmake')
    init_path = os.path.abspath(init_path)

    cwd = os.getcwd()
    for ii in range(pert_numb):
        path = os.path.join(out_dir, '01.pertubation')
        path = os.path.join(path, 'pert-%02d'%ii + '/tmp')
        print(' Working in {:s}'.format(path))
        create_path(path)
        shutil.copy2(incar, os.path.join(path, 'INCAR'))
        shutil.copy2(potcar, os.path.join(path, 'POTCAR'))
        #shutil.copy2(kpoints, os.path.join(path, 'KPOINTS'))
        os.chdir(path)
        for jj in range(neb_numb+2):
            neb_dir = '%02d'%jj
            poscar_dir = os.path.join(init_path, neb_dir)
            poscar = os.path.join(poscar_dir, 'POSCAR')
            work_path = os.path.join(path, neb_dir)
            create_path(work_path)
            if jj == 0 or jj == neb_numb+1:
                pos_in = poscar
                pos_out = os.path.join(work_path, 'POSCAR')
                shutil.copy2(pos_in, pos_out)
            else:
                pert_cmd = '/data/hutaiping/hutaiping/bin/newPy/geometry/tools/create_random_disturb.py'
                pert_cmd = 'python3 ' + pert_cmd + ' -etmax %f -ofmt vasp %s %d %f > /dev/null' %(pert_box, poscar, 1, pert_atom)
                sp.check_call(pert_cmd, shell=True)  #randomly perturbate for intermediate images
                pos_in = os.path.join(poscar_dir, 'POSCAR1.vasp')
                pos_out = os.path.join(work_path, 'POSCAR')
                shutil.copy2(pos_in, pos_out)
                os.remove(pos_in)
            os.chdir(cwd)
        
def structure_pertubation(json_file):
    jdata = json.load(open(json_file))

    pert_box = jdata['pert_box']
    pert_atom = jdata['pert_atom']
    pert_numb = jdata['pert_numb']

    incar = jdata['incar']
    potcar = jdata['potcar']
    #chgcar = jdata.get("chgcar")
    poscar = jdata["poscar"]

    cwd = os.getcwd()
    for ii in range(pert_numb):
        path = os.path.join("./", 'dpgen_work/NonSpinSP/struct_%02d'%ii)
        create_path(path)
        shutil.copy2(incar, os.path.join(path, 'INCAR'))
        shutil.copy2(potcar, os.path.join(path, 'POTCAR'))
        shutil.copy2(poscar, os.path.join(path, 'POSCAR'))
        #if chgcar is not None:
        #    shutil.copy2(incar, os.path.join(path, 'CHGCAR'))
        os.chdir(path)
        pert_cmd = '/root/DPAnalysis/tools/create_random_disturb.py'
        pert_cmd = 'python3 ' + pert_cmd + ' -etmax %f -ofmt vasp %s %d %f > /dev/null' %(pert_box, poscar, 1, pert_atom)
        sp.check_call(pert_cmd, shell=True)  #randomly perturbate for intermediate images
        #pos_in = os.path.abspath(os.path.join(path, 'POSCAR1.vasp'))
        #pos_out = os.path.abspath(os.path.join(path, 'POSCAR'))
        #shutil.copy2(pos_in, pos_out)
        shutil.copy2("POSCAR1.vasp", "POSCAR")
        os.remove("POSCAR1.vasp")
        os.chdir(cwd)


###########################################################


###########################################################
def dp_test_using_dp_calculator(iter_idx):
    import dpdata
    from ase.io import read
    import sys
    
    from deepmd.infer import DeepPot

    total_job_json = select_logs(make_iter_name(iter_idx)+'/02.fp', 'job.json')
    if len(total_job_json) == 0:
        raise ValueError(" Cannot find job.json file, please check again.")

    results = {}
    temps = []
    for job_json in total_job_json:
        temps.append(json.load(open(job_json))["temps"])
    temps_set = sorted(set(temps))
    for tmp in temps_set:
        results[tmp] = {"DFT": [], "DP":[]}

    for job_json in total_job_json:
        work_path = os.path.abspath(os.path.dirname(job_json) + os.path.sep + '.')
        outcar_file = os.path.join(work_path, 'OUTCAR')
        if os.path.exists(outcar_file):
            system = dpdata.LabeledSystem(outcar_file, fmt="vasp/outcar")
            atype = system["atom_types"]
            natoms = len(atype)
            results[json.load(open(job_json))["temps"]]["DFT"].append(system["energies"][0]/natoms)
            ase_read = read(outcar_file)
            coordinate = np.array(ase_read.get_positions()).reshape([1, -1])
            cell = np.array(ase_read.get_cell()).reshape([1, -1])
            dp = DeepPot("graph.000.pb")
            e, f, v = dp.eval(coordinate, cell, atype)
            results[json.load(open(job_json))["temps"]]["DP"].append(e[0][0]/natoms)
        else:
            print(" Cannot find the OUTCAR file in %s."%work_path)

    json_str = json.dumps(results, indent=4)
    with open("results_iter.%06d.json"%iter_idx, "w") as f:
        f.write(json_str)
###########################################################


###########################################################
def generate_anti_defect(ori_poscar, total_struct_numb, output_dir):
    """
    Generate anti defect structures.

    Written at 20211218
    """
    if not os.path.exists(ori_poscar):
        raise ValueError("The POSCAR file is not exist, please check again.") 

    import dpdata
    system = dpdata.System(ori_poscar, fmt="vasp/poscar")
    POTCAR_path = os.path.abspath("./POTCAR")
    INCAR_ISIF_3_path = os.path.abspath("./INCAR_ISIF_3")
    
    if not os.path.exists(POTCAR_path):
        raise ValueError("The POTCAR file is not exist, please check again.") 
    if not os.path.exists(INCAR_ISIF_3_path):
        raise ValueError("The INCAR file is not exist, please check again.")

    for counter in range(total_struct_numb):
        tmp_system = system.copy()
        tmp_system.replace("Li", "Co", 1)
        tmp_system.replace("Co", "Li", 1)
        output_path = os.path.join(output_dir, 'dpgen_work/struct_%d/0_ISIF_3'%counter)
        create_path(output_path)
        tmp_system.to_vasp_poscar(os.path.join(output_path, "POSCAR"))

        shutil.copy2(POTCAR_path, output_path)
        shutil.copy2(INCAR_ISIF_3_path, output_path+'/INCAR')
###########################################################


###########################################################
def generate_defect(ori_poscar, total_struc_numb, output_dir):
    import pymatgen
    POTCAR_path = os.path.abspath("./POTCAR")
    INCAR_ISIF_3_path = os.path.abspath("./INCAR_ISIF_3")

    if not os.path.exists(ori_poscar):
        raise ValueError("The POSCAR file is not exist, please check again.") 
    if not os.path.exists(POTCAR_path):
        raise ValueError("The POTCAR file is not exist, please check again.") 
    if not os.path.exists(INCAR_ISIF_3_path):
        raise ValueError("The INCAR file is not exist, please check again.")

    struct = pymatgen.Structure.from_file(ori_poscar)
    for counter in range(total_struc_numb):
        feed = random.sample(range(0, 12), 8)
        print(feed)
        feed_sort = arg_sort(feed, reverse=True)
        new_feed = [feed[idx] for idx in feed_sort] # sorted from large to small 
        
        tmp_struct = struct.copy()
        for idx in new_feed:
            tmp_struct.remove(tmp_struct[idx])
        output_path = os.path.join(output_dir, 'dpgen_work/struct_%d/0_ISIF_3'%counter)
        create_path(output_path)
        tmp_struct.to(filename=os.path.join(output_path, 'POSCAR'))

        shutil.copy2(POTCAR_path, output_path)
        shutil.copy2(INCAR_ISIF_3_path, output_path+'/INCAR')
###########################################################
