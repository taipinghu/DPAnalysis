#!/usr/bin/env python3

import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sb
import math
import scipy.constants as const
import json, glob, os
import subprocess as sp
from sklearn.metrics import mean_squared_error, r2_score

import warnings
warnings.filterwarnings('ignore')

import utils, fileio
from utils import file_to_list, make_iter_name, file_to_list
from fileio import parse_gmx_xvg, parse_deepmdkit_energy_out

from matplotlib.colors import LogNorm

colors=['red', 'blue', 'green']

###########################################################
##                    MD Analysis                        ##
###########################################################
def plot_RDF_comparison_figure(AIMD_xvg_files_list, DPMD_xvg_files_list, legend):
    fig = plt.figure()
    axes = fig.add_subplot()

    for itmp, aimd in enumerate(AIMD_xvg_files_list):
        x_data, y_data = parse_gmx_xvg(aimd, start_line=26, x_scale_factor=10)
        axes.plot(x_data, y_data, ls='--', color=colors[itmp])
    for itmp, dpmd in enumerate(DPMD_xvg_files_list):
        x_data, y_data = parse_gmx_xvg(dpmd, start_line=26, x_scale_factor=10)
        axes.plot(x_data, y_data, ls='-', color=colors[itmp])
    
    bold_axis_and_ticks(axes)
    set_legend(axes, label=legend, loc=1)
    return plt
###########################################################


###########################################################
def plot_five_figures_dp_vs_dft(energy_DFT: list, energy_DP: list, force_DFT: list, force_DP: list, prefix: str) -> None:
    """
    Plot five figures used in dp vs dft.
    Five Figures: Energy Comparison, Energy Error Distribution, Force Comparison, Force Error Distribution, Force Contour

    Args:
        energy_DFT: DFT-calculated energy
        energy_DP:  DP-calculated energy
        force_DFT:  DFT-calculated force
        force_DP:   DP-calculated force
        prefix: Figure name's prefix
    
    Returns:
        None, Five figures will be exported in current dir.

    Written at 20220119, by Taiping Hu.
    """
    plt = plot_energy_comparison(energy_DFT, energy_DP)
    plt.savefig("%s_Energy_Comparison.png"%prefix)
    plt.title("%s"%prefix)
    plt.close()
                        
    plt = plot_force_comparison(force_DFT, force_DP)
    plt.savefig("%s_Force_Comparison.png"%prefix)
    plt.title("%s"%prefix)
    plt.close()

    plt = plot_error_distribution(energy_DFT, energy_DP)
    plt.savefig("%s_Energy_Error_Distribution.png"%prefix)
    plt.title("%s"%prefix)
    plt.close()

    plt = plot_error_distribution(force_DFT, force_DP)
    plt.savefig("%s_Force_Error_Distribution.png"%prefix)
    plt.title("%s"%prefix)
    plt.close()

    plt = plot_force_contour_distribution(force_DFT, force_DP)
    plt.savefig("%s_Force_Error_Contour.png"%prefix) 
    plt.title("%s"%prefix)
    plt.close() 
    return None


def plot_energy_comparison(data_x, data_y, title=""):
    RMSE = np.sqrt(mean_squared_error(data_x, data_y))
    fig = plt.figure()
    axes = fig.add_subplot()
    axes.scatter(data_x, data_y, alpha=0.6)
    axes.plot((0, 1), (0, 1), transform=axes.transAxes,ls='--', c='b', label="1:1 line")
    bold_axis_and_ticks(axes)
    #set_legend(axes, label="RMSE:\n%.3e"%RMSE,loc=4)
    axes.set_aspect(1)
    axes.set_title("%s"%title, fontsize=20, fontweight="bold")

    return plt

def plot_force_comparison(data_x, data_y, title=""):
    fig = plt.figure()
    axes = fig.add_subplot()
    axes.scatter(data_x, data_y, alpha=0.6)
    axes.plot((0, 1), (0, 1), transform=axes.transAxes,ls='--', c='b', label="1:1 line")
    bold_axis_and_ticks(axes)
    axes.set_aspect(1)
    axes.set_title("%s"%title, fontsize=20, fontweight="bold")
    return plt

def plot_error_distribution(data_x, data_y, title=""):
    fig = plt.figure()
    axes = fig.add_subplot()
    energy_error = np.array(data_x) - np.array(data_y)
    sb.distplot(energy_error, hist=True, kde=True, bins=20, fit_kws=dict(linewidth=2.5))
    bold_axis_and_ticks(axes)
    axes.set_title("%s"%title, fontsize=20, fontweight="bold")
    return plt

def plot_force_contour_distribution(data_x, data_y, title=""):
    fig = plt.figure()
    axes = fig.add_subplot()
    total_x = np.array(data_x)
    total_y = np.array(data_y)            
    xmin, xmax = -6, 6
    ymin, ymax = -6, 6
    plt.hist2d(total_x, total_y, range=[[xmin, xmax], [ymin, ymax]],
                bins=400, norm=LogNorm(), cmap=plt.get_cmap('rainbow'))
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14, length=0)
    cb.ax.tick_params(which='minor', width=0)
    cb.outline.set_visible(False)
    plt.tick_params(width=2, labelsize=12)
    bold_axis_and_ticks(axes)
    axes.set_aspect(1)
    axes.set_title("%s"%title, fontsize=20, fontweight="bold")
    return plt
###########################################################


###########################################################
def plot_energy(file_name, atom_numb_list):
    data, pred = parse_deepmdkit_energy_out(file_name, atom_numb_list)
            
    fig = plt.figure()
    axes = fig.add_subplot()
    axes.set_aspect(1)
    axes.scatter(data, pred, alpha=0.6)
    axes.plot((0, 1), (0, 1), transform=axes.transAxes,ls='--', c='b', label="1:1 line")
    bold_axis_and_ticks(axes)

    return plt
###########################################################


###########################################################
def plot_energy_error_distribution(file_name, atom_numb_list):
    data, pred = parse_deepmdkit_energy_out(file_name, atom_numb_list)
    diff = data - pred
    fig = plt.figure()
    axes = fig.add_subplot()
    #axes.set_aspect(1)
    sb.distplot(diff, hist=True, kde=True, bins=20, fit_kws=dict(linewidth=2.5))
    bold_axis_and_ticks(axes)
#    set_legend(axes, ['',], loc=2)
    return plt
###########################################################


###########################################################
def plot_force_error_distribution(file_name):
    data_fx, pred_fx, data_fy, pred_fy, data_fz, pred_fz = fileio.data.parse_deepmdkit_force_out(file_name)
    diff_fx, diff_fy, diff_fz = data_fx - pred_fx, data_fy - pred_fy, data_fz - pred_fz
    fig = plt.figure()
    axes = fig.add_subplot()
    sb.distplot(diff_fx, hist=True, kde=True, bins=20, fit_kws=dict(linewidth=2.5))
    sb.distplot(diff_fy, hist=True, kde=True, bins=20, fit_kws=dict(linewidth=2.5))
    sb.distplot(diff_fz, hist=True, kde=True, bins=20, fit_kws=dict(linewidth=2.5))
    axes.set_aspect(1)
    bold_axis_and_ticks(axes)
    set_legend(axes, ['fx', 'fy', 'fz'], loc=1)
    return plt
###########################################################


###########################################################
def plot_force(file_name):
    data_fx, pred_fx, data_fy, pred_fy, data_fz, pred_fz = fileio.data.parse_deepmdkit_force_out(file_name)
    fig = plt.figure()
    xmin, xmax = -5, 5.1
    ymin, ymax = -5, 5.1

    axes = fig.add_subplot()
    axes.set_aspect(1)
    axes.scatter(data_fx, pred_fx, s=10, alpha=0.6)
    axes.scatter(data_fy, pred_fy, s=10, alpha=0.6)
    axes.scatter(data_fz, pred_fz, s=10, alpha=0.6)
    bold_axis_and_ticks(axes)
    set_legend(axes, label=['fx', 'fy', 'fz'], loc=2)
    axes.plot((0, 1), (0, 1), transform=axes.transAxes,ls='--', c='b', label="1:1 line")
    set_axis_limit_general(axes, xmin, xmax, ymin, ymax, x_step=5, y_step=5)
    
    return plt        
###########################################################


###########################################################
def plot_force_contour(file_name, bin=400):
    data_fx, pred_fx, data_fy, pred_fy, data_fz, pred_fz = fileio.data.parse_deepmdkit_force_out(file_name)

    total_x = np.array(data_fx.tolist() + data_fy.tolist() + data_fz.tolist())
    total_y = np.array(pred_fx.tolist() + pred_fy.tolist() + pred_fz.tolist())
    xmin = np.min(total_x)
    xmax = np.max(total_x)
    ymin = np.min(total_y)
    ymax = np.max(total_y)

    xmin, xmax = -6, 6
    ymin, ymax = -6, 6

    plt.hist2d(total_x, total_y, range=[[xmin, xmax], [ymin, ymax]],
                bins=bin, norm=LogNorm(), cmap=plt.get_cmap('rainbow'))
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14, length=0)
    cb.ax.tick_params(which='minor', width=0)
    cb.outline.set_visible(False)
    plt.tick_params(width=2, labelsize=12)

    ax = plt.gca()
    bold_axis_and_ticks(ax)

    return plt 
###########################################################


###########################################################
def plot_max_force_distribution_in_iteration(iter_index, param_file, template=False):
    """
    plot max force distribution in dpgen iteration.

    First written at 20210915
    Re-Written at 20211204
    """
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
###########################################################

def bold_axis_and_ticks(ax, x_label='', y_label='', linewidth=2, size=12):
    import matplotlib.ticker as mtick
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]
    ax.set_xlabel(x_label, fontsize=size, fontweight='bold')
    ax.set_ylabel(y_label, fontsize=size, fontweight='bold')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.2f"))
    return ax

def set_legend(ax, label, loc=8):
    font={'weight':'bold', 'size':12}
    ax.legend(label, prop=font, frameon=False, loc=loc)
    return ax

###########################################################
def set_axis_limit_general(ax, xmin, xmax, ymin, ymax, x_step, y_step):
    ax.xaxis.set_ticks(np.arange(float("{:.2f}".format(xmin)), float("{:.2f}".format(xmax)), x_step))       
    ax.yaxis.set_ticks(np.arange(float("{:.2f}".format(ymin)), float("{:.2f}".format(ymax)), y_step))
    ax.tick_params(labelsize=16)
    return ax
###########################################################

if __name__ == '__main__':
    #plt = plot_energy("model-000.e.out", [39, 39, 56, 48, 56, 48])
    plt = plot_max_force_distribution_in_iteration(7)  
    plt.show()
