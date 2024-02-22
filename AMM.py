# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:15:25 2024

@author: lafields2
"""

import pandas as pd
import csv

glycan_database_path = r"C:\Users\lawashburn\Documents\collaborations\Angel\glycan_database_formatted_long.csv"
glycan_database = pd.read_csv(glycan_database_path)

formatted_ms2_path = r"C:\Users\lawashburn\Documents\collaborations\Angel\PO_NP_heat_formatted.txt"

output_directory = r"C:\Users\lawashburn\Documents\collaborations\Angel"

z_states = [1,2,3]
H_mass = 1.00784
tolerance = 20

formatted_ms2 = pd.read_csv(formatted_ms2_path, sep=",",skiprows=[0], names= ['fragment_mz',
                                                                              'resolution',
                                                                              'fragment_charge',
                                                                              'fragment_intensity',
                                                                              'precursor_mz',
                                                                              'MS2_scan_number',
                                                                              'precursor_charge',
                                                                              'precursor_RT',
                                                                              'IonInjectTime',
                                                                              'ms1_scan',
                                                                              'precursor_intensity',
                                                                              'null'])


#%%
all_results_all_z = pd.DataFrame()

for z in z_states:
    formatted_ms2_filtered = formatted_ms2[formatted_ms2['precursor_charge'] == z]
    glycan_db_filtered = glycan_database[glycan_database['z'] == z]
    
    formatted_ms2_filtered['actual_monoisotopic_mass'] = (formatted_ms2_filtered['precursor_mz'] * formatted_ms2_filtered['precursor_charge']) - (H_mass * formatted_ms2_filtered['precursor_charge'])
    
    glycan_db_filtered['theoretical_monoisotopic_mass'] = (glycan_db_filtered['m/z'] * glycan_db_filtered['z']) - (H_mass * glycan_db_filtered['z'])
    
    tolerance_temp = tolerance *100
    
    formatted_ms2_filtered = formatted_ms2_filtered.sort_values(by='actual_monoisotopic_mass',ascending=True)
    glycan_db_filtered = glycan_db_filtered.sort_values(by='theoretical_monoisotopic_mass',ascending=True)
    
    merge_match2 = pd.merge_asof(formatted_ms2_filtered,glycan_db_filtered, left_on='actual_monoisotopic_mass', 
                                right_on='theoretical_monoisotopic_mass',
                                tolerance = tolerance_temp, allow_exact_matches=True,direction='forward') 
    
    merge_match3 = pd.merge_asof(formatted_ms2_filtered,glycan_db_filtered, left_on='actual_monoisotopic_mass', 
                                right_on='theoretical_monoisotopic_mass',
                                tolerance = tolerance_temp, allow_exact_matches=True,direction='backward') 
    
    all_matches = pd.concat([merge_match2,merge_match3])

    all_matches['Precursor error (ppm)'] = ((abs((all_matches['actual_monoisotopic_mass'])-
                                                          (all_matches['theoretical_monoisotopic_mass'])))/
                                                      (all_matches['actual_monoisotopic_mass'])) * 1E6
    
    all_matches_filtered = all_matches[all_matches['Precursor error (ppm)'] <= tolerance]
    all_matches_filtered = all_matches_filtered.drop(columns={'fragment_mz','resolution','fragment_intensity','fragment_charge','null'})
    all_matches_filtered_no_dups = all_matches_filtered.drop_duplicates()
    
    all_results_all_z = pd.concat([all_matches_filtered_no_dups,all_results_all_z])
    
output_path = output_directory + '\\precursor_AMM_results.csv'
with open(output_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        all_results_all_z.to_csv(filec,index=False)