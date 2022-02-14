#!/usr/bin/env python3

"""
Parse various file 
"""
import dpdata
import numpy as np
from utils import typemap_list_to_symbols, Molecule, file_to_list, process_bar
import os


###########################################################
##                  File I/O Functions                   ##
###########################################################
def read_file(file_or_dir_name, format):
    if format == "vasp/poscar":
        box, moles = read_vasp_poscar(file_or_dir_name)
    elif format == "lammps/dump":
        box, moles = read_lammps_multiframe_dump_file(file_or_dir_name)
    elif format == "vasp/contcar":
        box, moles = read_vasp_contcar(file_or_dir_name)
    return box, moles

def read_deepmd_npy(dir_name, format="deepmd/npy"):
    """
    read deepmd/npy file by calling dpdata
    """
    system = dpdata.LabeledSystem(dir_name, fmt=format)
    box, moles = {}, {}
    for iframe in range(system.get_nframes()):
        box[iframe], moles[iframe] = system["cells"][iframe], []
        natoms = len(system["coords"][iframe])
        atomic_symbols = typemap_list_to_symbols(system["atom_numbs"], system["atom_names"])
        for iatom in range(natoms):
            moles[iframe].append(Molecule([atomic_symbols[iatom],], [system["coords"][iframe][iatom],], [0.0, ]))
    return box, moles

def read_vasp_poscar(file_name, style="atom"):
    system = dpdata.System(file_name, fmt="vasp/poscar")

    box = system["cells"][0]
    moles = []
    atomic_symbols = typemap_list_to_symbols(system["atom_numbs"], system["atom_names"])

    return box, moles

def read_vasp_contcar(file_name, style="atom"):
    system = dpdata.System(file_name, fmt="vasp/poscar")

    box = system["cells"][0]
    moles = []
    atomic_symbols = typemap_list_to_symbols(system["atom_numbs"], system["atom_names"])

    return box, moles

def read_lammps_multiframe_dump_file(file_name):
    box, moles = {}, {}

    labeled_system = dpdata.System(file_name, fmt="lammps/dump")
    for iframe in range(labeled_system.get_nframes()):
        box = labeled_system["cells"][iframe]
        moles[iframe] = []

    return box, moles

def write_pdb(moles, file_name='out.pdb', mode=0, style='atoms'):
    """
    Write moles into pdb file.

    Args:
        moles       list        a list of molecule or Atoms objects.    
        mode        int         =0: output to a single file, =1: output to the separated files
    """

    if not isinstance(moles, (tuple, list)):
        raise TypeError('The argument moles must be iterable.')
    
    if style == 'atoms':  # in this case, the moles is a multi-frames Atoms object
        if mode == 0: # write into single file
            nframes = len(moles)
            with open(file_name, 'w') as f:
                for iframe in range(nframes):
                    f.write("MODEL         %r\n" %(iframe))     
                    f.write("TITLE     XXX t=  {:.5f}\n".format(iframe))
                    f.write("REMARK   Converted from other.write_pdb module\n")
                    cell = moles[iframe].get_cell_lengths_and_angles()
                    f.write('CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f}\n'.format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]))
                    natoms = len(moles[iframe].get_atomic_numbers())
                    symbols = moles[iframe].get_chemical_symbols()
                    coord = moles[iframe].get_positions()
                    for iatom in range(natoms):
                        f.write('%4s%7d%4s%5s%6d%4s%8.3f%8.3f%8.3f%6.2f%6.2f%12s\n'%("ATOM", iatom+1, symbols[iatom], "MOL", 1, '    ', coord[iatom][0], coord[iatom][1], coord[iatom][2], 1.0, 0.0, symbols[iatom]))
                    f.write('TER\n')
                    f.write('ENDMDL\n')
                    process_bar(iframe, nframes, " Write PDB file...")
        elif mode == 1: # write into separated file
            nframes = len(moles)
            for iframe in range(nframes):
                out_file = file_name.split('.')[0] + '-%d.pdb'%iframe
                with open(out_file, 'w') as f:
                    f.write("MODEL         %r\n" %(1))  
                    f.write("TITLE     XXX t=  {:.5f}\n".format(1)) 
                    f.write("REMARK   Converted from other.write_pdb module\n")
                    cell = moles[iframe].get_cell_lengths_and_angles()
                    f.write('CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f}\n'.format(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]))
                    natoms = len(moles[iframe].get_atomic_numbers())
                    symbols = moles[iframe].get_chemical_symbols()              
                    coord = moles[iframe].get_positions()
                    for iatom in range(natoms):
                        f.write('%4s%7d%4s%5s%6d%4s%8.3f%8.3f%8.3f%6.2f%6.2f%12s\n'%("ATOM", iatom+1, symbols[iatom], "MOL", 1, '    ', coord[iatom][0], coord[iatom][1], coord[iatom][2], 1.0, 0.0, symbols[iatom]))
                    f.write('TER\nENDMDL\n')    

    elif style == 'mole': # in this case, the moles is a single-frames Mole object      
        if mode == 0:
            nmoles = len(moles)
            with open(file_name, 'w') as f:
                f.write("MODEL         %r\n" %(1))      
                f.write("TITLE     XXX t=  {:.5f}\n".format(1))
                f.write("REMARK   Converted from other.write_pdb module\n")
                for imole in range(nmoles):
                    natoms = moles[imole].get_natoms()
                    symbols = moles[imole].get_atomic_symbols()
                    coord = moles[imole].get_atomic_coordinates()
                    for iatom in range(natoms):
                        f.write('%4s%7d%4s%5s%6d%4s%8.3f%8.3f%8.3f%6.2f%6.2f%12s\n'%("ATOM", iatom+1, symbols[iatom], "MOL", 1, '    ', coord[iatom][0], coord[iatom][1], coord[iatom][2], 1.0, 0.0, symbols[iatom]))
                f.write('TER\n')
                f.write('ENDMDL\n')
        elif mode == 1:
            pass
    else:
        pass
###########################################################


###########################################################
###                    deepmd-kit                       ###
###########################################################
def parse_deepmdkit_lcurve(file_name='lcurve.out'):
    """
    deepmd-kit version: 2.0.0, 
    #  step      rmse_trn    rmse_e_trn    rmse_f_trn         lr
    """
    lines = file_to_list(file_name)
    step = []
    rmse_e, rmse_f = [], []
    txt = lines[1:]
    for ii in range(0, len(txt), 1):
        tmp = txt[ii].split()
        step.append(int(tmp[0]))
        rmse_e.append(float(tmp[2]))
        rmse_f.append(float(tmp[3]))
    return step, rmse_e, rmse_f

def parse_deepmdkit_energy_out(file_name, atom_list):
    """
    Parse deepmdkit generated energy out file using -d.

    file_name       str     filename
    atom_list       list    atomic numbers

    Written at 20210908
    """
    lines = file_to_list(file_name)

    data_e, pred_e = [], []
    counter = -1
    for line in lines:
        tmp = line.split()
        if str(tmp[0]) == "#":
            counter += 1
            continue
        data_e.append(float(tmp[0])/atom_list[counter])
        pred_e.append(float(tmp[1])/atom_list[counter])
        
    data_e, pred_e = np.array(data_e), np.array(pred_e)
    return data_e, pred_e

def parse_deepmdkit_force_out(file_name):
    lines = file_to_list(file_name)
    data_fx, data_fy, data_fz = [], [], []
    pred_fx, pred_fy, pred_fz = [], [], []
    
    for line in lines[1:]:
        tmp = line.split()
        if str(tmp[0]) == "#":
            continue
        data_fx.append(float(tmp[0]))
        data_fy.append(float(tmp[1]))
        data_fz.append(float(tmp[2]))
        pred_fx.append(float(tmp[3]))
        pred_fy.append(float(tmp[4]))
        pred_fz.append(float(tmp[5]))
    data_fx, data_fy, data_fz = np.array(data_fx), np.array(data_fy), np.array(data_fz) 
    pred_fx, pred_fy, pred_fz = np.array(pred_fx), np.array(pred_fy), np.array(pred_fz) 

    return data_fx, pred_fx, data_fy, pred_fy, data_fz, pred_fz

def parse_model_devi_out(file_name='model_devi.out'):
    """
    Parse model devi output file.

    Written at 20210909
    """
    lines = file_to_list(file_name)
    
    max_devi_f = []
    for line in lines[1:]:
        tmp = line.split()
        max_devi_f.append(float(tmp[4]))
    
    return np.array(max_devi_f)

def parse_input_lammps(file_name='input.lammps'):

    lines = file_to_list(file_name)

    return float(lines[3].split()[3])

def parse_gmx_xvg(file_name, start_line=20, x_scale_factor=1.0, y_scale_factor=1.0):
    """
    parse gmx generated .xvg format file

    Args:
        start_line          int         data recorded start line
        scale_factor        float       unit convertion factor

    Written at 20210909
    """
    lines = file_to_list(file_name)

    x_data, y_data = [], []
    for line in lines[start_line:]:
        tmp = list(map(float, line.split()))
        x_data.append(tmp[0]*x_scale_factor)
        y_data.append(tmp[1]*y_scale_factor)
    
    x_data, y_data = np.array(x_data), np.array(y_data)
    return x_data, y_data
###########################################################    
