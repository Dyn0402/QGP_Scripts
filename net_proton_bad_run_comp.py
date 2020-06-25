#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 23 4:52 PM 2020
Created in PyCharm
Created as QGP_Scripts/net_proton_bad_run_comp.py

@author: Dylan Neff, dylan
"""

# No duplicates found but still remove just in case


def main():
    energies = [7, 11, 14, 19, 27, 39, 62]
    all_tosh_bad_runs = []
    all_note_bad_runs = []
    for energy in energies:
        print(f'\n{energy}GeV')
        tosh_bad_runs = get_tosh_bad_runs(energy)
        tosh_bad_dcaxy_runs = get_tosh_bad_dcaxy_runs(energy)
        note_bad_runs = get_note_bad_runs(energy)
        compare_bad_runs(tosh_bad_runs, note_bad_runs, tosh_bad_dcaxy_runs)

        if energy != 14:
            all_tosh_bad_runs += tosh_bad_runs
            all_note_bad_runs += note_bad_runs

    all_yang_bad_runs = get_yang_bad_runs()
    all_bad_ref_runs = get_ref_bad_runs()

    compare_all_bad_runs(all_tosh_bad_runs, all_note_bad_runs, all_yang_bad_runs, all_bad_ref_runs)

    print('donzo')


def get_tosh_bad_runs(energy):
    path = f'/home/dylan/Research/Results/6-23-20/RunNumber{energy}.h'
    bad_runs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        read_flag = False
        for line in lines:
            if "RunNumber[" in line:
                read_flag = True
                continue
            if read_flag:
                if '};' in line:
                    break
                line = line.strip().strip(',')
                bad_runs.append(int(line))
    bad_runs = list(dict.fromkeys(bad_runs))  # Remove duplicates

    return bad_runs


def get_tosh_bad_dcaxy_runs(energy):
    path = f'/home/dylan/Research/Results/6-23-20/BadEventBoundary{energy}.h'
    bad_runs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        read_flag = False
        for line in lines:
            if "BadRunFromDcaXY[" in line:
                read_flag = True
                continue
            if read_flag:
                if '};' in line:
                    break
                line = line.strip().strip(',').split(',')
                for run in line:
                    bad_runs.append(int(run))
    bad_runs = list(dict.fromkeys(bad_runs))  # Remove duplicates

    return bad_runs


def get_note_bad_runs(energy):
    path = '/home/dylan/Research/Results/6-23-20/analysis_note_bad_runs.txt'
    bad_runs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        read_flag = False
        for line in lines:
            if f'{{{energy},' in line:
                read_flag = True
                continue
            if read_flag:
                if '}}' in line:
                    break
                line = line.strip().strip('{').strip(',').split(',')
                for run in line:
                    bad_runs.append(int(run))
    bad_runs = list(dict.fromkeys(bad_runs))  # Remove duplicates

    return bad_runs


def get_yang_bad_runs():
    path = '/home/dylan/Research/Results/6-23-20/yangzz_bad_runs_3sig.txt'
    bad_runs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            bad_runs.append(int(line))
    bad_runs = list(dict.fromkeys(bad_runs))  # Remove duplicates

    return bad_runs


def get_ref_bad_runs():
    path = '/home/dylan/git/Research/QGP_Fluctuations/Tree_Reader/StRefMultCorr/BadRun.h'
    bad_runs = []
    with open(path, 'r') as file:
        lines = file.readlines()
        read_flag = False
        for line in lines:
            if 'const Int_t badrun_refmult' in line:
                if '2010' in line or '2011' in line:
                    read_flag = True
                    continue
            if read_flag:
                line = line.strip().strip('};').split(',')
                for run in line:
                    bad_runs.append(int(run))
                read_flag = False
    bad_runs = list(dict.fromkeys(bad_runs))  # Remove duplicates

    return bad_runs


def compare_bad_runs(tosh, note, dcaxy):
    tosh_not_note = [x for x in tosh if x not in note]
    note_not_tosh = [x for x in note if x not in tosh]
    dca_not_tosh = [x for x in dcaxy if x not in tosh]

    print(f'Length of Toshihiro: {len(tosh)}')
    print(f'Length of Analysis Note: {len(note)}')
    print(f'Length of DcaXY:: {len(dcaxy)}')
    print(f'Length of Toshihiro - Analysis Note: {len(tosh) - len(note)}')
    print(f'Length of Toshihiro - DcaXY: {len(tosh) - len(dcaxy)}')

    tosh_not_note_str = ''
    for run in tosh_not_note:
        tosh_not_note_str += str(run) + ','
    print(f'Toshihiro bad runs that are not in the analysis note {len(tosh_not_note)}: ')
    print(tosh_not_note_str)

    note_not_tosh_str = ''
    for run in note_not_tosh:
        note_not_tosh_str += str(run) + ','
    print(f'Analysis note bad runs that are not in Toshihiro\'s list {len(note_not_tosh)}: ')
    print(note_not_tosh_str)

    dca_not_tosh_str = ''
    for run in dca_not_tosh:
        dca_not_tosh_str += str(run) + ','
    print(f'DcaXY bad runs that are not in Toshihiro\'s list {len(dca_not_tosh)}: ')
    print(dca_not_tosh_str)


def compare_all_bad_runs(tosh, note, yang, ref):
    tosh = list(dict.fromkeys(tosh))  # Remove duplicates
    note = list(dict.fromkeys(note))  # Remove duplicates
    print('\n---------------------------------------------\n')
    print(f'Length of Toshihiro: {len(tosh)}')
    print(f'Length of Analysis Note: {len(note)}')
    print(f'Length of Yang: {len(yang)}')
    print(f'Length of Ref: {len(ref)}')

    tosh_not_note = [x for x in tosh if x not in note]
    note_not_tosh = [x for x in note if x not in tosh]
    yang_not_note = [x for x in yang if x not in note]
    note_not_yang = [x for x in note if x not in yang]
    ref_not_yang = [x for x in ref if x not in yang]
    yang_not_ref = [x for x in yang if x not in ref]

    yang_not_note_str = ''
    for run in yang_not_note:
        yang_not_note_str += str(run) + ','
    print(f'\nYang bad runs that are not in the analysis note {len(yang_not_note)}: ')
    print(yang_not_note_str)

    note_not_yang_str = ''
    for run in note_not_yang:
        note_not_yang_str += str(run) + ','
    print(f'\nAnalysis note bad runs that are not in Yang\'s list {len(note_not_yang)}: ')
    print(note_not_yang_str)

    yang_not_ref_str = ''
    for run in yang_not_ref:
        yang_not_ref_str += str(run) + ','
    print(f'\nYang bad runs that are not in ref {len(yang_not_ref)}: ')
    print(yang_not_ref_str)

    ref_not_yang_str = ''
    for run in ref_not_yang:
        ref_not_yang_str += str(run) + ','
    print(f'\nRef bad runs that are not in Yang\'s list {len(ref_not_yang)}: ')
    print(ref_not_yang_str)


if __name__ == '__main__':
    main()
