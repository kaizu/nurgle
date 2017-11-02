#!/usr/bin/env python
# -*- coding: utf-8 -*-

from logging import basicConfig, getLogger, DEBUG
log_ = getLogger(__name__)
import os
import os.path
import json
import csv
import functools
import re

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

def assume_source(data, default_location='CCO-CYTOSOL'):
    for i in range(len(data['reactions'])):
        reaction = data['reactions'][i]
        if reaction['id'].startswith('__') or reaction['id'] == 'MetaFlux_obj':
            log_.debug('reaction [{}] has no corresponding reaction.'.format(reaction['id']))
            continue
        elif 'name' not in reaction:
            log_.warn('reaction [{}] has no name.'.format(reaction['id']))
            continue

        if '/' not in reaction['name']:
            rid = reaction['name']
        else:
            rid = reaction['name'].split('/')[0]
            for name in reaction['metabolites']:
                if re.sub(r'[-\+]', '_', rid).endswith('_{}'.format(name[: -2])):
                    rid = rid[: -len(name)+1]
                    break
            else:
                log_.warn('reaction [{}] has reactants [{}].'.format(reaction['name'].split('/')[0].replace('-', '_'), [name for name, coef in reaction['metabolites'].items()]))
                continue

        loc = default_location
        mobj = re.search(r'\[([^\[\]]+)\]', rid)
        if mobj is not None:
            loc = mobj[1]
            rid = rid[: -len(loc)-2]
        assert '[' not in rid and ']' not in rid

        reaction['source'] = dict(id=rid, location=loc)
        log_.debug('reaction id of [{}] was assumed to be [{}@{}].'.format(reaction['id'], rid, loc))

    return data

def generate_ecocyc_fba(ECOCYC_VERSION="21.1", showall=False):
    log_.info('generate_ecocyc_fba(ECOCYC_VERSION="{}")'.format(ECOCYC_VERSION))

    filename = os.path.join(INPUTS_PATH, ECOCYC_VERSION, 'data/fba/fba-examples/ecocyc-21.0-gem-cs-glucose-tea-none.json')
    if not os.path.isfile(filename):
        raise RuntimeError('an input file [{}] could not be found'.format(filename))
    log_.info('read a file [{}]'.format(filename))
    with open(filename, 'r') as fin:
        data = json.load(fin)

    data = assume_source(data)
    data = solve_fba_using_cobra(data)

    compounds = []
    for reaction in data['reactions']:
        if not showall and reaction['flux'] == 0.0:
            continue
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
            is_constant = 0
            writer.writerow((name, val, is_constant))

    filename = os.path.join(OUTPUTS_PATH, 'metabolism.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')
        for reaction in data['reactions']:
            if not showall and reaction['flux'] == 0.0:
                continue
            assert ';' not in reaction['id'] and ':' not in reaction['id']
            assert all(';' not in name and ':' not in name for name in reaction['metabolites'])
            flux = reaction['flux']
            # reversible = (reaction.get('lower_bound', 0.0) < 0.0)

            reactants = ';'.join('{}:{}'.format(name, -coef) for name, coef in reaction['metabolites'].items() if coef < 0)
            products = ';'.join('{}:{}'.format(name, coef) for name, coef in reaction['metabolites'].items() if coef > 0)

            # vf = functools.reduce(lambda x, y: (x * (compounds[y[0]] ** -y[1])) if y[1] < 0 else x, reaction['metabolites'].items(), 1.0)
            # vr = functools.reduce(lambda x, y: (x * (compounds[y[0]] ** y[1])) if y[1] > 0 else x, reaction['metabolites'].items(), 1.0)
            vf = functools.reduce(lambda x, y: (x * compounds[y[0]]) if y[1] < 0 else x, reaction['metabolites'].items(), 1.0)
            vr = functools.reduce(lambda x, y: (x * compounds[y[0]]) if y[1] > 0 else x, reaction['metabolites'].items(), 1.0)
            vfmax, vrmax = 0.0, 0.0

            if flux > 0.0:
                # assert vf > 0.0
                # vfmax = flux / vf
                assert vf > 0.0 and vr > 0.0
                vfmax = 2 * flux / vf
                vrmax = flux / vr
            elif flux < 0.0:
                # assert vr > 0.0
                # vrmax = -flux / vr
                assert vf > 0.0 and vr > 0.0
                vfmax = -flux / vf
                vrmax = 2 * -flux / vr
            else:
                # flux == 0.0
                if reaction['lower_bound'] >= 0 or reaction['upper_bound'] <= 0:
                    pass  # do nothing
                else:
                    vfmax = 1.0 / vf
                    vrmax = 1.0 / vr

            writer.writerow((reaction['id'], reactants, products, vfmax, vrmax))


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
