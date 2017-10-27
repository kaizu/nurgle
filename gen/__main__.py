#!/usr/bin/env python
# -*- coding: utf-8 -*-

from logging import basicConfig, getLogger, DEBUG
log_ = getLogger(__name__)
import os
import os.path
import json
import csv
import functools

INPUTS_PATH = 'inputs'
OUTPUTS_PATH = 'outputs'


def solve_fba_using_cobra(data):
    log_.info('solve_fba_using_cobra')

    try:
        import cobra
    except ImportError as e:
        log_.error('module [cobra] is required.')
        raise e

    log_.info('cobra.__version__ = {}'.format(cobra.__version__))

    cmodel = cobra.Model(data['name'])

    for reaction in data['reactions']:
        r = cobra.Reaction(reaction['id'])
        r.add_metabolites(dict((cobra.Metabolite(name), coef) for name, coef in reaction['metabolites'].items()))
        r.lower_bound = reaction.get('lower_bound', 0)
        r.upper_bound = reaction.get('upper_bound', 3000)
        cmodel.add_reaction(r)

    cmodel.objective = "__SEC_EX_obj_met_c"
    soln = cmodel.optimize()

    for i in range(len(data['reactions'])):
        reaction = data['reactions'][i]
        reaction['flux'] = soln.x_dict[reaction['id']]

    return data

def generate_ecocyc_fba(ECOCYC_VERSION="21.1"):
    log_.info('generate_ecocyc_fba(ECOCYC_VERSION="{}")'.format(ECOCYC_VERSION))

    filename = os.path.join(INPUTS_PATH, ECOCYC_VERSION, 'data/fba/fba-examples/ecocyc-21.0-gem-cs-glucose-tea-none.json')
    if not os.path.isfile(filename):
        raise RuntimeError('an input file [{}] could not be found'.format(filename))
    log_.info('read a file [{}]'.format(filename))
    with open(filename, 'r') as fin:
        data = json.load(fin)

    data = solve_fba_using_cobra(data)

    compounds = []
    for reaction in data['reactions']:
        compounds.extend(reaction['metabolites'])
    default = 1.0
    compounds = dict((name, default) for name in set(compounds))

    filename = os.path.join(INPUTS_PATH, 'compounds.csv')
    if os.path.isfile(filename):
        log_.info('read a file [{}]'.format(filename))
        with open(filename, 'r') as fin:
            for i, row in enumerate(csv.reader(fin, lineterminator='\n')):
                if len(row) != 2:
                    raise RuntimeError('invalid format at Line {}. the size of a row must be 2. {} was given'.format(i + 1, len(row)))
                compounds[row[0]] = float(row[1])

    filename = os.path.join(OUTPUTS_PATH, 'compounds.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')
        for name, val in compounds.items():
            writer.writerow((name, val))

    filename = os.path.join(OUTPUTS_PATH, 'metabolism.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')
        for reaction in data['reactions']:
            assert ';' not in reaction['id'] and ':' not in reaction['id']
            assert all(';' not in name and ':' not in name for name in reaction['metabolites'])
            metabolites = ';'.join('{}:{}'.format(name, coef) for name, coef in reaction['metabolites'].items())
            denom = functools.reduce(lambda x, y: x * (compounds[y[0]] ** y[1]), reaction['metabolites'].items(), 1.0)
            if denom != 0.0:
                writer.writerow((reaction['id'], metabolites, reaction['flux'] / denom))
            else:
                if reaction['flux'] == 0.0:
                    log_.warning('a flux of reaction [{}] is 0. vmax is 1.0 as its default'.format(reaction['id']))
                else:
                    log_.error('a flux of reaction [{}] is 0, but the one at equilibrium is non-zero [{}]'.format(reaction['id'], reaction['flux']))
                writer.writerow((reaction['id'], metabolites, 1.0))


if __name__ == "__main__":
    basicConfig(level=DEBUG)

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='a path to input files', default='inputs')
    parser.add_argument('-o', '--output', type=str, help='a path to output files', default='outputs')
    args = parser.parse_args()

    INPUTS_PATH = args.input
    OUTPUTS_PATH = args.output

    if not os.path.isdir(OUTPUTS_PATH):
        log_.info('make a directory [{}]'.format(OUTPUTS_PATH))
        os.mkdir(OUTPUTS_PATH)

    generate_ecocyc_fba()
