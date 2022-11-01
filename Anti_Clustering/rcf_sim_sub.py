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

from calc_binom_slices import find_sim_sets


def main():
    submit_xml_path = '/star/u/dneff/git/QGP_Fluctuations/Tree_Reader/subs/submit_sub.xml'
    sets = get_sets()
    for set_i in sets:
        submit_set(set_i, submit_xml_path)
    print('donzo')


def get_sets():
    clust_alg = 'anticlmulti'
    names = ['set_group', 'set_name']
    existing_sets_path = '/star/u/dneff/gpfs/tree_reader_data/Data_Sim/'
    existing_sets = find_sim_sets(existing_sets_path, ['flat80', clust_alg, 'resample'], ['test'], False)
    existing_sets = None if existing_sets.size == 0 else existing_sets
    print(f'Existing Sets:\n{existing_sets}')

    sets = [
        # ['flat80_anticlmulti_spread4_amp05_resample', 'Sim_spread4_amp05_flat80_anticlmulti_norotate_resample_'],
        # ['flat80_anticlmulti_spread4_amp06_resample', 'Sim_spread4_amp06_flat80_anticlmulti_norotate_resample_']
    ]

    # Run specific (amp, spread) pairs
    # sa_pairs = [('15', '35'), ('15', '0225'), ('15', '1'), ('1', '35'), ('1', '005'), ('1', '125'), ('225', '045'),
    #             ('225', '4'), ('225', '08'), ('225', '35'), ('225', '05'), ('225', '5'), ('225', '035'),
    #             ('225', '0225'), ('225', '005'), ('225', '1'), ('225', '175'), ('25', '0175'), ('25', '4'),
    #             ('25', '08'), ('25', '35'), ('25', '05'), ('25', '5'), ('25', '035'), ('25', '0225'), ('25', '005'),
    #             ('25', '1'), ('25', '175'), ('275', '045'), ('275', '0175'), ('275', '35'), ('275', '5'),
    #             ('275', '035'), ('275', '0225'), ('275', '005'), ('275', '1'), ('275', '175'), ('2', '0175'),
    #             ('2', '08'), ('2', '35'), ('2', '035'), ('2', '005'), ('2', '1'), ('2', '175'), ('325', '045'),
    #             ('325', '0'), ('325', '01'), ('325', '5'), ('325', '035'), ('325', '0225'), ('325', '005'),
    #             ('325', '07'), ('325', '1'), ('325', '175'), ('35', '045'), ('35', '0'), ('35', '015'), ('35', '05'),
    #             ('35', '5'), ('35', '0225'), ('35', '005'), ('35', '1'), ('35', '175'), ('375', '045'),
    #             ('375', '0175'), ('375', '15'), ('375', '0'), ('375', '015'), ('375', '5'), ('375', '035'),
    #             ('375', '225'), ('375', '0225'), ('375', '07'), ('375', '1'), ('375', '175'), ('375', '3'),
    #             ('3', '045'), ('3', '35'), ('3', '05'), ('3', '5'), ('3', '035'), ('3', '0225'), ('3', '005'),
    #             ('3', '07'), ('3', '1'), ('3', '175'), ('3', '3'), ('45', '045'), ('45', '09'), ('45', '4'),
    #             ('45', '15'), ('45', '25'), ('45', '35'), ('45', '0'), ('45', '015'), ('45', '05'), ('45', '2'),
    #             ('45', '5'), ('45', '06'), ('45', '225'), ('45', '0225'), ('45', '07'), ('45', '1'), ('4', '045'),
    #             ('4', '09'), ('4', '03'), ('4', '125'), ('4', '25'), ('4', '35'), ('4', '0'), ('4', '015'), ('4', '5'),
    #             ('4', '2'), ('4', '035'), ('4', '225'), ('4', '0225'), ('4', '07'), ('4', '1'), ('4', '175'),
    #             ('4', '3'), ('5', '045'), ('5', '09'), ('5', '0175'), ('5', '03'), ('5', '25'), ('5', '0'),
    #             ('5', '01'), ('5', '015'), ('5', '035'), ('5', '2'), ('5', '225'), ('5', '06'), ('5', '0225'),
    #             ('5', '07'), ('5', '1'), ('5', '3')]
    # for spread, amp in sa_pairs:
    #     sets.append([f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
    #                  f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_'])

    # Rerun all spreads/amps
    amps = ['0', '002', '004', '005', '006', '008', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04',
            '045', '05', '06', '07', '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5',
            '6', '7', '8', '9', '99']
    spreads = ['001', '01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25',
               '275', '3', '325', '35', '375', '4']  # '45', '5']  # Run 45 and 5 on Old PC
    for amp in amps:
        for spread in spreads:
            if existing_sets is None or \
                    existing_sets[(existing_sets['spread'] == spread) & (existing_sets['amp'] == amp)].size == 0:
                sets.append([f'flat80_{clust_alg}_spread{spread}_amp{amp}_resample',
                             f'Sim_spread{spread}_amp{amp}_flat80_{clust_alg}_norotate_resample_'])

    # # New spreads all amps
    # amps = ['0', '005', '01', '0125', '015', '0175', '02', '0225', '025', '03', '035', '04', '045', '05', '06', '07',
    #         '08', '09', '1', '125', '15', '175', '2', '225', '25', '3', '35', '4', '45', '5', '6', '7', '8', '9', '99']
    # spreads = ['06', '065', '07', '075', '08', '085', '09']
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_{clust_alg}_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_{clust_alg}_norotate_resample_'])
    #
    # # New amps all spreads
    # amps = ['0125', '6', '7', '8', '9', '99']
    # spreads = ['001', '01', '02', '05', '06', '065', '07', '075', '08', '085', '09', '1', '15', '2', '225', '25',
    #            '275', '3', '325', '35', '375', '4']  # '45', '5']  # Run 45 and 5 on Old PC
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_{clust_alg}_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_{clust_alg}_norotate_resample_'])

    # # New amps missed spreads
    # amps = ['275', '3', '325', '35', '375', '4', '45', '5']
    # spreads = ['001', '01']
    # for amp in amps:
    #     for spread in spreads:
    #         sets.append([f'flat80_anticlmulti_spread{spread}_amp{amp}_resample',
    #                      f'Sim_spread{spread}_amp{amp}_flat80_anticlmulti_norotate_resample_'])

    # for seti in sets:
    #     print(seti)
    print(f'Number of sets to resubmit: {len(sets)}')

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

    sleep(1)

    os.remove(new_xml_path)


if __name__ == '__main__':
    main()
