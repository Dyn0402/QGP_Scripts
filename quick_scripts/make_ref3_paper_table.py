#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 04 1:23 PM 2025
Created in PyCharm
Created as QGP_Scripts/make_ref3_paper_table.py

@author: Dylan Neff, Dylan
"""

import re
import pandas as pd


def main():
    """
    Read refmult and refmult3 bin edges from Param.h for each energy and create a table for the ref3 paper.
    :return:
    """
    strefcor_param_path = 'C:/Users/Dylan/source/repos/Dyn0402/QGP_Fluctuations/Tree_Reader/StRefMultCorr/Param.h'
    df_ref3 = extract_strefcorr_param_data(strefcor_param_path)
    include_cols = ['Energy', 'Year', 'Bin Edges']
    exclude_energies = ['200', '54', '14.5']
    df_ref3 = df_ref3[~df_ref3['Energy'].isin(exclude_energies)]
    # Sort on integer of energy then integer of trigger start
    df_ref3['Energy_float'] = df_ref3['Energy'].astype(float)
    df_ref3['Trigger Start Int'] = df_ref3['Trigger Start'].astype(int)
    df_ref3 = df_ref3.sort_values(['Energy_float', 'Trigger Start Int'])
    # For each energy choose only the lowest trigger start
    df_ref3 = df_ref3.groupby('Energy_float').head(1)
    write_latex_table(df_ref3, include_cols=include_cols)
    print('donzo')


def write_latex_table(df, file_path='ref3_table.txt', include_cols=None):
    """
    Write ref3 bin edges to a latex table. Include columns for energy, year, trigger start, trigger end,
    vz start, vz end, and bin edges.
    :param df:
    :param file_path:
    :param include_cols:
    :return:
    """
    if include_cols is None:
        include_cols = ['Energy', 'Year', 'Trigger Start', 'Trigger End', 'Vz Start', 'Vz End', 'Bin Edges']
    with open('ref3_table.txt', 'w') as file:
        file.write('\\begin{table*}[]\n')
        file.write('\\centering\n')
        file.write('\\begin{tabular}{|c|c|c|c|c|c|c|}\n')
        file.write('\\hline\n')
        file.write(' & '.join(include_cols) + ' \\\\\n')
        file.write('\\hline\n')
        for i, row in df.iterrows():
            file.write(' & '.join([str(row[col]) for col in include_cols]) + ' \\\\\n')
        file.write('\\hline\n')
        file.write('\\end{tabular}\n')
        file.write('\\end{table*}')


def extract_strefcorr_param_data(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    # Extract content between "mParamStr_ref3" and "};"
    match = re.search(r'mParamStr_ref3(.*?)\};', content, re.DOTALL)
    if not match:
        print("Error: mParamStr_ref3 section not found in the file.")
        pd.DataFrame()

    section_content = match.group(1).strip()
    section_content = section_content[section_content.index('\n')+1:]  # Remove first line
    # section_content = section_content[section_content.index('\n')+1:]  # Remove second line
    section_content = section_content[:section_content.rindex('\n')]  # Remove last line
    section_content = section_content.replace('\t', '')
    section_content = section_content.replace('"', '')
    section_content = section_content.replace(' ', '')

    entries = section_content.split('},\n{')

    extracted_data = []
    for entry in entries:
        lines = entry.split('\n')  # Split by newline markers
        if len(lines) >= 2:
            year, energy, trigs, vz = lines[1].strip(',').split(':')
            trig_start, trig_end = trigs.split(',')
            vz_start, vz_end = vz.split(',')
            bin_edges = lines[2].strip(',').split(',')
            bin_edges = [int(edge) for edge in bin_edges]
            bin_edges_str = ', '.join([str(edge) for edge in bin_edges])
            extracted_data.append([year, energy, trig_start, trig_end, vz_start, vz_end, bin_edges_str])

    # Convert to Pandas DataFrame
    df = pd.DataFrame(extracted_data, columns=['Year', 'Energy', 'Trigger Start', 'Trigger End', 'Vz Start', 'Vz End', 'Bin Edges'])
    return df

if __name__ == '__main__':
    main()
