#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on December 27 2:02 PM 2021
Created in PyCharm
Created as QGP_Scripts/rcf_sim_sub

@author: Dylan Neff, Dyn04
"""

import os
from time import sleep


def main():
    submit_xml_path = '/star/u/dneff/git/QGP_Fluctuations/Tree_Reader/subs/submit_sub.xml'
    sets = get_sets()
    for set_i in sets:
        submit_set(set_i, submit_xml_path)
    print('donzo')


def get_sets():
    names = ['set_group', 'set_name']
    sets = [
        # ['flat80_anticlmulti_spread4_amp05_resample', 'Sim_spread4_amp05_flat80_anticlmulti_norotate_resample_'],
        # ['flat80_anticlmulti_spread4_amp06_resample', 'Sim_spread4_amp06_flat80_anticlmulti_norotate_resample_']
    ]

    # Run specific (amp, spread) pairs


    # Rerun all spreads/amps
    # amps = ['0', '005', '01', '015', '0175', '02', '0225', '025', '03', '035', '04', '045', '05', '06', '07', '08',
    #         '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5']
    # spreads = ['001', '01', '02', '05', '1', '15', '2', '225', '25', '275', '3', '325', '35', '375', '4', '45', '5']
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_'])

    # New spreads all amps
    # amps = ['0', '005', '01', '015', '02', '025', '03', '035', '04', '045', '05', '06', '07', '08', '09',
    #         '12', '15', '175', '2', '225', '25']
    # spreads = ['001', '01', '02', '05', '1', '15', '2', '225', '25', '275', '3', '325', '35', '375', '4', '45', '5']
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_'])

    # # New amps all spreads
    # amps = ['275', '3', '325', '35', '375', '4', '45', '5']
    # spreads = ['001', '01', '02', '05', '1', '15', '2', '225', '25', '275', '3', '325', '35', '375', '4', '45', '5']
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_'])

    # # New amps missed spreads
    # amps = ['275', '3', '325', '35', '375', '4', '45', '5']
    # spreads = ['001', '01']
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_'])

    # for seti in sets:
    #     print(seti)
    # print(len(sets))

    sets = [dict(zip(names, x)) for x in sets]

    return sets


def submit_set(set_i, xml_path):
    with open(xml_path, 'r') as file:
        lines = file.readlines()

    new_lines = []
    for line in lines:
        if './Release/Tree_Reader' in line:
            new_lines.append(f'\t\t./Release/Tree_Reader {set_i["set_group"]} {set_i["set_name"]}')
        else:
            new_lines.append(line)

    new_xml_path = f'sub_{set_i["set_name"]}.xml'
    with open(new_xml_path, 'w') as file:
        file.writelines(new_lines)

    print(f'star-submit {new_xml_path}')
    os.system(f'star-submit {new_xml_path}')

    sleep(2)

    os.remove(new_xml_path)


if __name__ == '__main__':
    main()
