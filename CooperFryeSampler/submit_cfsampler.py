#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 22 4:30 PM 2023
Created in PyCharm
Created as QGP_Scripts/submit_cfsampler

@author: Dylan Neff, Dylan
"""

import os


def main():
    n_jobs = 20
    submit = False
    energies = [7, 19, 27, 39, 62]
    run_types = {'CooperFrye_all': 'input.AuAu.7.7.C0-5',
                 'CooperFrye_EV_all': 'input.AuAu.7.7.C0-5.EVHRG',
                 'CooperFrye_EVb342_all': 'input.AuAu.7.7.C0-5.EVHRG_b342'}

    for energy in energies:
        for output_dir, input_template in run_types.items():
            replacements = {'CooperFrye_protons': output_dir,
                            'nProcesses="200"': f'nProcesses="{n_jobs}"',
                            'input.AuAu.7.7.C0-5.EVHRG': {input_template}}
            template_path = f'/star/u/dneff/git/CooperFryeSamplerRunner/subs/sub{energy}.xml'
            out_path = f'/star/u/dneff/gpfs/data/CooperFrye/gen_subs/sub_{output_dir}_{energy}GeV.xml'

            gen_xml(template_path, out_path, replacements)

            if submit:
                os.system(f'star-submit {out_path}')

    print('donzo')


def gen_xml(template_path, out_path, replacements):
    with open(template_path, 'r') as in_file:
        text = in_file.read()
    for orig, new in replacements.items():
        text = text.replace(orig, new)
    with open(out_path, 'w') as out_file:
        out_file.write(text)


if __name__ == '__main__':
    main()
