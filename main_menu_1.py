#!/usr/bin/env python3

import shutil
import pandas as pd
import os, sys
import numpy as np
import glob
import json

import warnings
warnings.filterwarnings('ignore')

import dpdata
from re import template
from utils import string_2_array, figure_utils, exported_figure_file, figure_output_options
from utils import select_logs, make_iter_name, file_to_list
from utils import select_dirs, select_all_dirs, process_bar, create_path, process_bar

from utils import DPAlog

train = "00.train"
model_devi = '01.model_devi'
fp = "02.fp"

###########################################################
def dpgen_iteration_functions():
    while True:
        print("\n")
        print("========================================================")
        print("                    Main Menu 1                        ")
        print("                                                       ")        
        print(" (  2)  Print three ratios in dpgen iteration          ")
        print(" (  3)  Calculate trust level in dpgen iteration       ")
        print(" (  4)  Calculate trust level in init data             ")
        print("                                                       ")
        print(" (  9)  Calculate the cost of init bulk                ")
        print(" ( 10)  Calculate the cost of dpgen iteration          ")
        print("                                                       ")
        print("                Model Deviation Analysis               ")
        print(" ( 11)  Plot max force error distribution with different temperaturs")
        print(" ( 12)  Plot max force error distribution with different pressures")
        print(" ( 13)  Plot model devi distribution for single out file")
        print("                                                       ")
        print("               Data Collection Functions               ")
        print(" ( 90)  Collect different systems from init data       ")
        print(" ( 91)  Collect different systems from iteration data  ")
        print(" ( 92)  Collect different scale from init data         ")
        print(" ( 93)  Collect different initial structure data from init")
        print(" ( 94)  Collect different data with different temperatures from iteration data")
        print(" ( 95)  Collect all data from init and iteration       ")
        print(" ( 96)  Collect all data in current dir from OUTCAR file")
        
        print("                                                       ")
        print(" (100)  Remove traj files in dpgen iteration           ")
        print(" (101)  Remove redundant ckpt files in train           ")
        print(" (102)  Remove redundant files generated by dpdispatcher")
        print("\n Tips: Input -10 to return main menu                 ")
        print("========================================================")
        jsel = input("\n Please input the menu index: ")

        if jsel == str(11) or jsel == str(12):
            str_input = input("\n If use template? y(Y) or n(N)\n")
            if_template = False
            if str_input == "y" or str_input == "Y":
                if_template = True
            str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
            param_file = input("\n Please input param file used by dpgen.\n")
            counter = 1
            for iter_idx in string_2_array(str_in):
                if jsel == str(11):
                    plt = plot_max_force_distribution_with_different_temperatures(iter_idx, param_file=param_file, template=if_template)
                    plt.savefig("Temps_Max_Force_Error_Distribution_Iter_%06d.png"%iter_idx)
                elif jsel == str(12):
                    plt = plot_max_force_distribution_with_different_press(iter_idx, param_file=param_file, template=if_template)
                    plt.savefig("Press_Max_Force_Error_Distribution_Iter_%06d.png"%iter_idx)
                #plt.show()
                process_bar(counter, len(string_2_array(str_in)), "Plot figures...")
                counter += 1
            plt.close()
        elif jsel == str(13):
            model_devi_file = input("\n Please input the model devi file name: ")
            plt = plot_model_devi_distribution_for_single_file(model_devi_file)
            plt.show()
        elif jsel == str(2):
            ksel = int(input("\n Print (1) total ratio or (2) individual ratio?\n"))
            if ksel == 1:
                str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
                print(" Iteration Index     Candidate          Failed        Accurate")
                DPAlog.info(" Iteration Index     Candidate          Failed        Accurate")
                for iter_idx in string_2_array(str_in):
                    candidate, accurate, failed = print_three_total_ratios_in_dpgen_iteration(iter_idx)
                    print("    {:6d}       {:10.2f} %     {:10.2f} %    {:10.2f} %".format(iter_idx, candidate*100, failed*100, accurate*100))
                    DPAlog.info("    {:6d}       {:10.2f} %     {:10.2f} %    {:10.2f} %".format(iter_idx, candidate*100, failed*100, accurate*100))
            elif ksel == 2:
                str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
                for iter_idx in string_2_array(str_in):
                    print_three_ratios_in_dpgen_iteration(iter_idx)
            print("\n Hints: Those results has been exported into DPAanalysis.log file.")
        elif jsel == str(3):
            print("\n Notes: All the results will be exported to DPAnalysis.log file. ")
            str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
            for iter_idx in string_2_array(str_in):
                calculate_trust_level_in_dpgen_iteration(iter_idx)
        elif jsel == str(4):
            str_in = input("\n Please input the dir root name include init data: ")
            calculate_trust_level_in_init_data(str_in)
        elif jsel == str(9):
            print(" Hints: First, you should put all the init tasks in current dir, then the program will automatically calculate all init tasks.")
            print(" By default, the cost will be exported into the time_init.xlsx file.")
            calculate_cost_init("./")
        elif jsel == str(10):
            print(" Hints: You should run this script in your dpgen iteration dir.")
            print(" By default, the cost will be exported into the time_iter.xlsx file.")
            calculate_cost_iteration("./")
        elif jsel == str(-10):
            break
        elif jsel == str(100):
            str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
            for iter_idx in string_2_array(str_in):
                counter = 0
                traj_path = glob.glob(os.path.join(make_iter_name(iter_idx), "01.model_devi/task.*/traj"))
                traj_path.sort()
                if len(traj_path) == 0:
                    print(" Cannot find traj dirs in iter.%06d..."%iter_idx)
                else:
                    for traj in traj_path:
                        shutil.rmtree(traj)
                        process_bar(counter, len(traj_path), " Remove traj file in iter.%06d..."%iter_idx)
                        counter += 1
        elif jsel == str(90):
            print("\n Extract different systems from init data. ")
            str_in = input("\n Please input the root dir name of init data: ")
            extract_different_system_from_init_data(str_in)
        elif jsel == str(91):
            print("\n Extract different systems from iteration data.")
            str_in = input("\n Please input the root dir name of init data: ")
            extract_different_system_from_iteration_data(str_in)
        elif jsel == str(92):
            str_in = input("\n Please input the init param json file: ")
            extract_different_scale_from_init_data(json_file=str_in)
        elif jsel == str(93):
            extract_different_AIMD_traj_from_init_data()
        elif jsel == str(94):
            collect_different_data_with_different_temps_from_iteration("./")
        elif jsel == str(95):
            str_in = input("\n Please input the param_run json file name, default is param_run.json (press Enter): ")
            if str_in == "":
                collect_all_data_from_init_and_iter_data(json_file="param_run.json")       
            else:
                collect_all_data_from_init_and_iter_data(json_file=str_in)
        elif jsel == str(96):
            root_dir = input("\n Please input the root dir name: ")
            collect_all_data_in_current_dir_from_all_OUTCAR_file(root_dir=root_dir)
        elif jsel == str(101):
            str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
            stop_batch = int(input("\n Please input the stop batch number: "))
            for iter_idx in string_2_array(str_in):
                train_path_total = glob.glob(os.path.join(make_iter_name(iter_idx), "00.train/00*"))
                train_path_total.sort()
                for train_path in train_path_total:
                    for numb in range(stop_batch-4000, stop_batch, 1000):
                        ckpt_files = glob.glob(os.path.join(train_path, "model.ckpt-%d.*"%numb))
                        ckpt_files.sort()
                        for ckpt_file in ckpt_files:
                            if os.path.exists(ckpt_file):
                                os.remove(ckpt_file)
        elif jsel == str(102):
            str_in = input("\n Please input the iter index, e.g., 1, 3-6, 8, 10-11, denotes 1, 3, 4, 5, 6, 8, 10, 11.\n")
            for iter_idx in string_2_array(str_in):
                counter = 1
                redundant_file = []
                iter_path = make_iter_name(iter_idx)
                redundant_file.extend(glob.glob(os.path.join(iter_path, "00.train/*.sub")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "01.model_devi/*.sub")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "02.fp/*.sub")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "00.train/*fail")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "01.model_devi/*fail")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "02.fp/*fail")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "00.train/*tag_finished")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "01.model_devi/*tag_finished")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "02.fp/*tag_finished")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "00.train/backup")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "01.model_devi/backup")))
                redundant_file.extend(glob.glob(os.path.join(iter_path, "02.fp/backup")))
                if len(redundant_file) == 0:
                    print(" Cannot find the file in %s..."%iter_path)
                else:
                    for file_tmp in redundant_file:
                        if os.path.isdir(file_tmp):
                            shutil.rmtree(file_tmp)
                        else:
                            os.remove(file_tmp)
                        process_bar(counter, len(redundant_file), " Delete file in %s..."%iter_path)
        else:
            print(" !!! Unspported menu index, please input again. !!!")
###########################################################


###########################################################
##              functions of data collection             ##
###########################################################
def collect_all_data_from_init_and_iter_data(json_file: str) -> None:
    import dpdata

    jdata = json.load(open(json_file))

    total_data = {}
    total_system = []
    init_data_prefix = jdata.get('init_data_prefix', '')
    init_data_sys = jdata.get('init_data_sys', [])
    for tmp_init in init_data_sys:
        labeled_system = dpdata.LabeledSystem(os.path.join(init_data_prefix, tmp_init), fmt="deepmd/npy")
        total_system.append(labeled_system.formula)

    model_devi_jobs = jdata.get('model_devi_jobs', {})
    numb_jobs = len(model_devi_jobs)
    for iter_idx in range(numb_jobs):
        iter_data = glob.glob(os.path.join(make_iter_name(iter_idx), '02.fp', 'data.[0-9]*[0-9]'))
        iter_data.sort()
        for tmp_iter in iter_data:
            labeled_system = dpdata.LabeledSystem(tmp_iter, fmt = 'deepmd/npy')
            total_system.append(labeled_system.formula)
    total_system = set(total_system)
    
    for system in total_system:
        total_data[system] = None

    for tmp_init in init_data_sys:
        labeled_system = dpdata.LabeledSystem(os.path.join(init_data_prefix, tmp_init), fmt="deepmd/npy")
        for system in total_system:
            if labeled_system.formula == system:
                if total_data[system] is None:
                    total_data[system] = labeled_system
                else:
                    total_data[system].append(labeled_system)
    for iter_idx in range(numb_jobs):
        iter_data = glob.glob(os.path.join(make_iter_name(iter_idx), '02.fp', 'data.[0-9]*[0-9]'))
        iter_data.sort()
        for tmp_iter in iter_data:
            labeled_system = dpdata.LabeledSystem(tmp_iter, fmt = 'deepmd/npy')
            for system in total_system:
                if labeled_system.formula == system:
                    if total_data[system] is None:
                        total_data[system] = labeled_system
                    else:
                        total_data[system].append(labeled_system)
    str_in = input("\n Please input the dir name to export: ")
    for system in total_system:
        total_data[system].to_deepmd_raw("%s/%s"%(str_in, system))
        total_data[system].to_deepmd_npy("%s/%s"%(str_in, system))


def collect_all_data_in_current_dir_from_all_OUTCAR_file(root_dir):
    import dpdata
    all_outcar = select_logs(root_dir, "OUTCAR")

    final_sys = None

    for outcar in all_outcar:
        labeled_system = dpdata.LabeledSystem(outcar, fmt="vasp/outcar")
        if final_sys is None:
            final_sys = labeled_system[-1]
        else:
            final_sys.append(labeled_system[-1])
    final_sys.to_deepmd_raw("deepmd")
    final_sys.to_deepmd_npy("deepmd")

def get_all_temps_in_all_iteration(root_dir) -> list:
    temps = []
 
    total_job_json = glob.glob(os.path.join(root_dir, "iter.*/02.fp/task.*/job.json"))
    for job_json in total_job_json:
        temps.append(json.load(open(job_json))["temps"])
    return sorted(set(temps))

def collect_different_data_with_different_temps_from_iteration(root_dir):
    from dpdata import MultiSystems
    fp_path_total = glob.glob(os.path.join(root_dir, "iter.*/02.fp/task.*"))
    fp_path_total.sort()

    total_temps = get_all_temps_in_all_iteration(root_dir)
    result = {}
    for temp in total_temps:
        result[temp] = MultiSystems()
    counter = 1
    for fp_path in fp_path_total:
        job_json = os.path.join(fp_path, "job.json")
        outcar_path = os.path.join(fp_path, "OUTCAR")
        if not os.path.exists(job_json):
            #print(" job.json file is not exist in %s..."%fp_path)
            continue
        if not os.path.exists(outcar_path):
            #print(" OUTCAR file is not exist in %s..."%fp_path)
            continue
        jdata = json.load(open(job_json))
        labeled_system = dpdata.LabeledSystem(outcar_path, fmt="vasp/outcar")
        for temp in total_temps:
            if temp in jdata.values():
                if result[temp] is None:
                    result[temp] = labeled_system
                else:
                    result[temp].append(labeled_system)
        process_bar(counter, len(fp_path_total))
        counter += 1
    print(result)

def extract_different_AIMD_traj_from_init_data(scale=1.000):
    import dpdata

    fp_path_total = glob.glob(os.path.join("./", "POSCAR.*/02.md/sys-*/*/*/OUTCAR"))
    fp_path_total.sort()

    create_path("Init_Data")
    counter = 0
    for fp_path in fp_path_total:
        labeled_system = dpdata.LabeledSystem(fp_path, fmt="vasp/outcar")
        labeled_system.to_deepmd_raw("Init_Data/traj.%03d"%counter)
        labeled_system.to_deepmd_npy("Init_Data/traj.%03d"%counter, set_size=len(labeled_system))
        counter += 1
        process_bar(counter, len(fp_path_total))

def extract_different_scale_from_init_data(json_file: str, infomode=0) -> None:
    import dpdata

    jdata = json.load(open(json_file))
    scale_factor = jdata["scale"]
    if infomode == 0:
        print(" Scale factor: {}".format(scale_factor))
    
    result = {}
    for scale in scale_factor:
        fp_path_total = glob.glob(os.path.join("./", "POSCAR.*/02.md/sys-*/scale-%.3f/*/OUTCAR"%scale))
        fp_path_total.sort()
        result[scale] = None
        for fp_path in fp_path_total:
            labeled_system = dpdata.LabeledSystem(fp_path, fmt="vasp/outcar")
            if result[scale] is None:
                result[scale] = labeled_system
            else:
                result[scale].append(labeled_system)
        result[scale].to_deepmd_raw("Init_Data/scale-%.3f"%scale)
        result[scale].to_deepmd_npy("Init_Data/scale-%.3f"%scale, set_size=len(result[scale]))

def extract_different_system_from_init_data(init_dir:str, infomode=0) -> None:
    import dpdata

    total_dirs = select_all_dirs(init_dir, "deepmd")

    total_system = []
    for dir_name in total_dirs:
        labeled_system = dpdata.LabeledSystem(dir_name, fmt="deepmd/npy")
        total_system.append(labeled_system.formula)
    total_system = set(total_system)
    
    
    result = {}
    for tmp_sys in total_system:
        result[tmp_sys] = None
    
    for dir_name in total_dirs:
        labeled_system = dpdata.LabeledSystem(dir_name, fmt="deepmd/npy")
        for tmp_sys in total_system:
            if labeled_system.formula == tmp_sys:
                if result[tmp_sys] is None:
                    result[tmp_sys] = labeled_system
                else:
                    result[tmp_sys].append(labeled_system)
    for tmp_sys in total_system:
        result[tmp_sys].to_deepmd_raw("Init_Data/%s"%tmp_sys)
        result[tmp_sys].to_deepmd_npy("Init_Data/%s"%tmp_sys, set_size=len(result[tmp_sys]))
        if infomode == 0:
            print(" Exporting %s, number of frames: %d"%(tmp_sys, len(result[tmp_sys]))) 

def extract_different_system_from_iteration_data(root_dir, infomode=0):
    import dpdata

    fp_path_total = glob.glob(os.path.join(root_dir, "iter.*/02.fp/data.*"))
    fp_path_total.sort()

    system_total = []
    counter = 1
    for fp_path in fp_path_total:
        labeled_system = dpdata.LabeledSystem(fp_path, fmt="deepmd/npy")
        system_total.append(labeled_system.formula)
        process_bar(counter, len(fp_path_total))
        counter += 1
    system_total = set(system_total)

    print("\n Now, all data has been loaded, exporting those data based on different system.\n")
    result = {}
    for sys_tmp in system_total:
        result[sys_tmp] = None
   
    for fp_path in fp_path_total:
        labeled_system = dpdata.LabeledSystem(fp_path, fmt="deepmd/npy")
        for tmp_sys in system_total:
            if labeled_system.formula == tmp_sys:
                if result[tmp_sys] is None:
                    result[tmp_sys] = labeled_system
                else:
                    result[tmp_sys].append(labeled_system)
    for tmp_sys in system_total:
        result[tmp_sys].to_deepmd_raw("Iter_Data/%s"%tmp_sys)
        result[tmp_sys].to_deepmd_npy("Iter_Data/%s"%tmp_sys, set_size=len(result[tmp_sys]))
        if infomode == 0:
            print(" Exporting %s, number of frames: %d"%(tmp_sys, len(result[tmp_sys])))
###########################################################


###########################################################
##        Time cost scripts are adopted from dfz         ##
###########################################################
def calculate_cost_iteration(root_iter):
    """
    This code is written by dfz

    Adopted at 20211217
    """     
    iters = [os.path.join(root_iter, p) 
             for p in os.listdir(root_iter) if 'iter' in p]
    time_iter = pd.DataFrame(index = list(range(len(iters)*3)),
                             columns = ['iter', 'task', 'num',
                                        'maxt', 'mint', 'meant',
                                        'price_per_unit', 'price'])
    
    task_info = {'00.train':'train.log',
                 '01.model_devi':'model_devi.log',
                 '02.fp':'OUTCAR'}
    
    i = -1
    for piter in iters:
        print(piter[-6::])
        for ptask in task_info.keys():
            i += 1
            time_iter.loc[i, 'iter'] = int(piter[-6::])
            time_iter.loc[i, 'task'] = ptask
            local_root = os.path.join(piter, ptask)
            if not os.path.exists(local_root):
                continue
            print(local_root)
            logs = select_logs(local_root, fname = task_info[ptask])
            time_iter.loc[i, 'num'] = len(logs)
            
            if int(ptask[:2:]) == 0:
                times  = [get_time_train(log) for log in logs]
            elif int(ptask[:2:]) == 1:
                times  = [get_time_devi(log) for log in logs]
            else:
                times  = [get_time_fp(log) for log in logs]
                
            times = np.array(times, dtype=np.float64)
            times = times[~np.isnan(times)]
            if times.size == 0:
                continue
            time_iter.loc[i, 'maxt'] = times.max()
            time_iter.loc[i, 'mint'] = times.min()
            time_iter.loc[i, 'meant'] = times.mean()

    time_iter.sort_values(['task', 'iter'], inplace=True)
    time_iter.to_excel('time_iter.xlsx')

    return None

def calculate_cost_init(root_init):
    """
    This code is written by dfz

    Adopted at 20211217
    """
    logs = select_logs(root_init, fname='OUTCAR')
    time_init = pd.DataFrame(index = list(range(len(logs))),
                             columns = ['filename', 'time in h',
                                        'price_per_unit', 'price'])
    time_init['filename'] = logs
    t = [get_time_fp(f) for f in logs]
    time_init['time in h'] = t
    time_init.to_excel('time_init.xlsx')
    return None

def get_time_devi(filename):
    """
    This code is written by dfz

    Adopted at 20211217
    """
    with open(filename, "r") as fip:
        lines = fip.readlines()
    for ii in lines[::-1]:
        if "Total wall time" in ii:
            tmp = ii.split()[3].split(":")
            return float(tmp[0]) + float(tmp[1])/60 + float(tmp[2])/3600

def get_time_train(filename):
    """
    This code is written by dfz

    Adopted at 20211217
    """    
    with open(filename, "r") as fip:
        lines = fip.readlines()
    for ii in lines[::-1]:
        if "wall time" in ii:
            return float(ii.split()[4])/3600

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
def print_three_total_ratios_in_dpgen_iteration(iter_index):
    """
    Print three total ratios in dpgen iteration.

    Written by Taiping Hu at 20211211.
    """
    from utils import get_file_line_number
    fp_dir = os.path.join(make_iter_name(iter_index), fp)
    candidate_out_total = glob.glob(os.path.join(fp_dir, "candidate.*"))
    accurate_out_total = glob.glob(os.path.join(fp_dir, "rest_accurate.*"))
    failed_out_total = glob.glob(os.path.join(fp_dir, "rest_failed.*"))
    candidate_out_total.sort()
    accurate_out_total.sort()
    failed_out_total.sort()
    
    candidate_numb, accurate_numb, failed_numb = [], [], []
    for candidate_out in candidate_out_total:
        candidate_numb.append(get_file_line_number(candidate_out))
    for accurate_out in accurate_out_total:
        accurate_numb.append(get_file_line_number(accurate_out))
    for failed_out in failed_out_total:
        failed_numb.append(get_file_line_number(failed_out))

    candidate, accurate, failed = np.sum(candidate_numb), np.sum(accurate_numb), np.sum(failed_numb)
    total = candidate + accurate + failed
    candidate_ratio, accurate_ratio, failed_ratio = candidate / total, accurate / total, failed / total
    return candidate_ratio, accurate_ratio, failed_ratio
###########################################################


###########################################################
def print_three_ratios_in_dpgen_iteration(iter_index):
    """
    Print three ratios for each system in dpgen iteration.

    Written by Taiping Hu at 20211211.
    """
    from utils import get_file_line_number
    fp_dir = os.path.join(make_iter_name(iter_index), fp)
    candidate_out_total = glob.glob(os.path.join(fp_dir, "candidate.*"))
    candidate_out_total.sort()

    system_index = []
    for candidate_out in candidate_out_total:
        system_index.append(candidate_out.split('.')[-2])

    print(" Iteration Index: %d"%iter_index)
    for sys in system_index:
        candidate_out = os.path.join(fp_dir, 'candidate.shuffled.%s.out'%sys)
        accurate_out  = os.path.join(fp_dir, 'rest_accurate.shuffled.%s.out'%sys)
        failed_out    = os.path.join(fp_dir, 'rest_failed.shuffled.%s.out'%sys)
        candidate_numb = get_file_line_number(candidate_out)
        accurate_numb  = get_file_line_number(accurate_out)
        failed_numb    = get_file_line_number(failed_out)
        total = candidate_numb + accurate_numb + failed_numb
        print(" System {0:s} {1:9s} : {2:10d} in {3:10d} {4:6.2f} %".format(sys, "candidate", candidate_numb, total, candidate_numb/total*100))
        print(" System {0:s} {1:9s} : {2:10d} in {3:10d} {4:6.2f} %".format(sys, "failed", failed_numb, total, failed_numb/total*100))
        print(" System {0:s} {1:9s} : {2:10d} in {3:10d} {4:6.2f} %".format(sys, "accurate", accurate_numb, total, accurate_numb/total*100))
        DPAlog.info(" System {0:s} {1:9s} : {2:10d} in {3:10d} {4:6.2f} %".format(sys, "candidate", candidate_numb, total, candidate_numb/total*100))
        DPAlog.info(" System {0:s} {1:9s} : {2:10d} in {3:10d} {4:6.2f} %".format(sys, "failed", failed_numb, total, failed_numb/total*100))
        DPAlog.info(" System {0:s} {1:9s} : {2:10d} in {3:10d} {4:6.2f} %".format(sys, "accurate", accurate_numb, total, accurate_numb/total*100))
###########################################################


###########################################################
def calculate_trust_level_in_dpgen_iteration(iter_idx):
    fp_task_path = glob.glob(os.path.join(make_iter_name(iter_idx), "02.fp/task.*"))
    fp_task_path.sort()

    system_index = []
    for fp_path in fp_task_path:
        system_index.append(os.path.basename(fp_path).split(".")[1])
    system_index.sort()
    set_tmp = set(system_index)
    system_index = list(set_tmp)
    system_index.sort()

    for sys in system_index:
        data_paths = glob.glob(os.path.join(make_iter_name(iter_idx), "02.fp/data.%s"%sys))
        data_paths.sort()
        for data_path in data_paths:
            force_raw_file = os.path.join(data_path, "force.raw")
            if not os.path.exists(force_raw_file):
                print(" Cannot find the force.raw file in %s"%data_path)
                continue
            lines = file_to_list(force_raw_file)
            data = []
            for line in lines:
                tmp = np.abs(list(map(float, line.split())))
                data.extend(tmp)
            DPAlog.info(" %s/02.fp/data.%s : %.4f  %.4f"%(make_iter_name(iter_idx), sys, 0.2*np.mean(data), 0.5*np.mean(data)))
            print(" %s/02.fp/data.%s : %.4f  %.4f"%(make_iter_name(iter_idx), sys, 0.2*np.mean(data), 0.5*np.mean(data)))
###########################################################


###########################################################
def calculate_trust_level_in_init_data(dir_name):
    total_force_raw_files = sorted(select_logs(dir_name, "force.raw"))

    for force_raw_file in total_force_raw_files:
        lines = file_to_list(force_raw_file)
        data = []
        for line in lines:
            tmp = np.abs(list(map(float, line.split())))
            data.extend(tmp)
        print("%s: %.4f  %.4f"%(force_raw_file, 0.2*np.mean(data), 0.5*np.mean(data)))
###########################################################


###########################################################
##              plot model devi distribution             ##
###########################################################
def plot_max_force_distribution_with_different_temperatures(iter_index, param_file, template=False):
    import matplotlib.pyplot as plt
    import seaborn as sb
    from plot_new import bold_axis_and_ticks, set_legend
    jdata = json.load(open(param_file))

    model_devi_jobs = jdata['model_devi_jobs']
    if template == True:
        temps = model_devi_jobs[iter_index]["rev_mat"]["lmp"]["V_TEMP"]
        press = model_devi_jobs[iter_index]["rev_mat"]["lmp"]["V_PRES"]
    else:
        temps = model_devi_jobs[iter_index]['temps']
        press = model_devi_jobs[iter_index]['press']
    trust_lo = jdata['model_devi_f_trust_lo']
    trust_hi = jdata['model_devi_f_trust_hi']

    counter = 0
    temps_dir = {}
    for temp in temps:
        temps_dir[temp] = []

    work_path = os.path.join(make_iter_name(iter_index), '01.model_devi')
    work_path_glob = glob.glob(os.path.join(os.path.abspath(work_path), "task.*.*"))
    work_path_glob.sort()
    for work_tmp in work_path_glob:
        for temp in temps:
            lines = file_to_list(work_tmp+'/input.lammps')
            if template:
                tmp_str = "variable        TEMP            equal %d"%temp
            else:
                tmp_str = "variable        TEMP            equal %f"%temp
            if tmp_str in lines:
                temps_dir[temp].append(work_tmp+'/model_devi.out')
    
    model_devi_out_data = {}
    for tmp_temp in temps:
        model_devi_out_data[tmp_temp] = []
        for tmp_file in temps_dir[tmp_temp]:
            all_conf = np.loadtxt(tmp_file)
            for ii in range(all_conf.shape[0]):
                model_devi_out_data[tmp_temp].append(all_conf[ii][4])

    fig = plt.figure()
    axes = fig.add_subplot()
    for tmp_temp in temps:
        sb.distplot(model_devi_out_data[tmp_temp], hist=None, kde=True, bins=20, fit_kws=dict(linewidth=2.5))   
        plt.xlim(xmin=0, xmax=0.25)
        plt.ylim(ymin=0, ymax=100)

    plt.axvline(x=trust_lo, c='black', ls=':')  
    plt.axvline(x=trust_hi, c='black', ls=':')
    
    bold_axis_and_ticks(axes)
    set_legend(axes, [str(tmp)+'K' for tmp in model_devi_out_data.keys()], loc=1)
    plt.title('Iteration %d'%iter_index, fontweight="bold", fontsize=20)
    return plt

def plot_max_force_distribution_with_different_press(iter_index, param_file, template=False):
    import matplotlib.pyplot as plt
    import seaborn as sb
    from plot_new import bold_axis_and_ticks, set_legend
    jdata = json.load(open(param_file))

    model_devi_jobs = jdata['model_devi_jobs']
    if template == True:
        temps = model_devi_jobs[iter_index]["rev_mat"]["lmp"]["V_TEMP"]
        press = model_devi_jobs[iter_index]["rev_mat"]["lmp"]["V_PRES"]
    else:
        temps = model_devi_jobs[iter_index]['temps']
        press = model_devi_jobs[iter_index]['press']
    trust_lo = jdata['model_devi_f_trust_lo']
    trust_hi = jdata['model_devi_f_trust_hi']

    counter = 0
    press_dir = {}
    for pres in press:
        press_dir[pres] = []

    work_path = os.path.join(make_iter_name(iter_index), '01.model_devi')
    work_path_glob = glob.glob(os.path.join(os.path.abspath(work_path), "task.*.*"))
    work_path_glob.sort()
    for work_tmp in work_path_glob:
        for pres in press:
            lines = file_to_list(work_tmp+'/input.lammps')
            if template:
                tmp_str = "variable        PRES            equal %d"%pres
            else:
                tmp_str = "variable        PRES            equal %f"%pres
            if tmp_str in lines:
                press_dir[pres].append(work_tmp+'/model_devi.out')
    
    model_devi_out_data = {}
    for tmp_pres in press:
        model_devi_out_data[tmp_pres] = []
        for tmp_file in press_dir[tmp_pres]:
            all_conf = np.loadtxt(tmp_file)
            for ii in range(all_conf.shape[0]):
                model_devi_out_data[tmp_pres].append(all_conf[ii][4])

    fig = plt.figure()
    axes = fig.add_subplot()
    for tmp_pres in press:
        sb.distplot(model_devi_out_data[tmp_pres], hist=None, kde=True, bins=20, fit_kws=dict(linewidth=2.5))   
        plt.xlim(xmin=0, xmax=0.25)
        plt.ylim(ymin=0, ymax=100)

    plt.axvline(x=trust_lo, c='black', ls=':')  
    plt.axvline(x=trust_hi, c='black', ls=':')

    bold_axis_and_ticks(axes)
    set_legend(axes, [str(tmp)+'Bar' for tmp in model_devi_out_data.keys()], loc=1)
    plt.title('Iteration %d'%iter_index, fontweight="bold", fontsize=20)
    return plt

def plot_model_devi_distribution_for_single_file(model_devi_file):
    import matplotlib.pyplot as plt
    import seaborn as sb
    from plot_new import bold_axis_and_ticks, set_legend
    model_devi_out_data = []
    all_conf = np.loadtxt(model_devi_file)
    for ii in range(all_conf.shape[0]):
        model_devi_out_data.append(all_conf[ii][4])
    
    fig = plt.figure()
    axes = fig.add_subplot()
    
    sb.distplot(model_devi_out_data, hist=None, kde=True, bins=20, fit_kws=dict(linewidth=2.5))   
    #plt.xlim(xmin=0, xmax=0.25)
    #plt.ylim(ymin=0, ymax=100)
    bold_axis_and_ticks(axes)
    return plt
###########################################################