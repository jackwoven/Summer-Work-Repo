import pymatgen
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester

def effective_mass_val(array):
    deltax = 1 
    first_d = -(array[:-2]-array[2:])
    second_d = (array[:-2]-array[1:-1]+array[2:])
    
    effective_masses=[]
    
    for i in range(len(first_d) -1):
        if np.sign(first_d[i]) == 1:
            if np.sign(first_d[i+1]) == -1:
                m = (first_d[i+1] - first_d[i])/deltax
                point = i/deltax - first_d[i]/m
                
                effective_masses.append(i)
        elif np.sign(first_d[i]) == 0:
            if np.sign(first_d[i+1]) == -1:
                effective_masses.append(i)
    return effective_masses

def effective_mass_cond(array):
    deltax = 1 
    first_d = -(array[:-2]-array[2:])
    second_d = (array[:-2]-array[1:-1]+array[2:])
    
    effective_masses=[]
    
    for i in range(len(first_d) -1):
        if np.sign(first_d[i]) == -1:
            if np.sign(first_d[i+1]) == 1:
                m = (first_d[i+1] - first_d[i])/deltax
                point = i/deltax - first_d[i]/m
                
                effective_masses.append(i)
        elif np.sign(first_d[i]) == 0:
            if np.sign(first_d[i+1]) == 1:
                effective_masses.append(i)
    return effective_masses
    

def min_finder(material_id):
    
    with MPRester() as mpr:
        bs = mpr.get_bandstructure_by_material_id(material_id)

    for i in bs.bands.items():
        aa =i[1] 

    cond_band_array = []
    val_band_array = []

    for j in aa:
        if  bs.efermi + 5.5 > min(j) > bs.efermi:
            cond_band_array.append(j)
        elif  bs.efermi > min(j) > bs.efermi-5.5:
            val_band_array.append(j)

    new_cond_band_array = np.array(cond_band_array)
    new_val_band_array = np.array(val_band_array)
   
    bs_dict = bs.as_dict()
    
    val_band = new_val_band_array
    cond_band = new_cond_band_array

    val_list1 = effective_mass_val(val_band[0])
    val_list2 = effective_mass_val(val_band[1])
    
    val_band_mins_1 = []
    val_band_mins_2 = []
    
    for i in val_list1:
        index = i
        energy = val_band[0][i]
        emass_left = 1/(val_band[0][i-2]-2*val_band[0][i-1]+val_band[0][i])
        emass_right = 1/(val_band[0][i+2]-2*val_band[0][i+1]+val_band[0][i])
        kpoint = bs_dict['kpoints'][index]
        val_band_mins_1.append([index, energy, emass_left, emass_right, kpoint])
    
    for j in val_list2:
        index = j
        energy = val_band[1][j]
        emass_left = 1/(val_band[1][j-2]-2*val_band[1][j-1]+val_band[1][j])
        emass_right = 1/(val_band[1][j+2]-2*val_band[1][j+1]+val_band[1][j])
        kpoint = bs_dict['kpoints'][index]
        val_band_mins_2.append([index, energy, emass_left, emass_right, kpoint])
    
    
    cond_list1 = effective_mass_cond(cond_band[0])
    cond_list2 = effective_mass_cond(cond_band[1])
    
    cond_band_mins_1 = []
    cond_band_mins_2 = []
    
    for k in cond_list1:
        index = k
        energy = cond_band[0][k]
        emass_left = 1/(cond_band[0][k-2]-2*cond_band[0][k-1]+cond_band[0][k])
        emass_right = 1/(cond_band[0][k+2]-2*cond_band[0][k+1]+cond_band[0][k])
        kpoint = bs_dict['kpoints'][index]
        cond_band_mins_1.append([index, energy, emass_left, emass_right, kpoint])
    
    for m in cond_list2:
        index = m
        energy = cond_band[1][m]
        emass_left = 1/(cond_band[1][m-2]-2*cond_band[1][m-1]+cond_band[1][m])
        emass_right = 1/(cond_band[1][m+2]-2*cond_band[1][m+1]+cond_band[1][m])
        kpoint = bs_dict['kpoints'][index]
        cond_band_mins_2.append([index, energy, emass_left, emass_right, kpoint])
    
    return val_band_mins_1, val_band_mins_2, cond_band_mins_1, cond_band_mins_2


min_finder("mp-567636")