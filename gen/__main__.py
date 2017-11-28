#!/usr/bin/env python
# -*- coding: utf-8 -*-

from logging import basicConfig, getLogger, DEBUG, INFO
log_ = getLogger(__name__)
import os
import os.path
import json
import csv
import functools
import re
from collections import defaultdict
from itertools import combinations
import enum
from math import isclose, floor
import copy

import numpy

INPUTS_PATH = 'inputs'
OUTPUTS_PATH = 'outputs'

REL_TOL = 1e-7
ABS_TOL = 1e-12

class Direction(enum.Enum):
    FORWARD = enum.auto()
    REVERSE = enum.auto()
    REVERSIBLE = enum.auto()
    DAMMED = enum.auto()


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

    reaction_ids = dict((reaction['id'], i) for i, reaction in enumerate(data['reactions']))

    for i in range(len(data['reactions'])):
        reaction = data['reactions'][i]
        reaction['flux'] = soln.x_dict[reaction['id']]

        lower = reaction.get('lower_bound', 0)
        upper = reaction.get('upper_bound', 3000)
        assert upper >= lower
        if upper == lower:
            reaction['direction'] = Direction.DAMMED
            assert isclose(reaction['flux'], lower, rel_tol=REL_TOL, abs_tol=ABS_TOL)
        elif lower >= 0:
            reaction['direction'] = Direction.FORWARD
            if reaction['flux'] < 0.0: print(reaction['flux'])
            assert reaction['flux'] > -ABS_TOL
        elif upper <= 0:
            if reaction['flux'] > 0.0: print(reaction['flux'])
            reaction['direction'] = Direction.REVERSE
            assert reaction['flux'] < +ABS_TOL
        else:
            reaction['direction'] = Direction.REVERSIBLE

        # if reaction['id'].startswith('__NUT_EX_'):
        #     if reaction['direction'] == Direction.DAMMED:
        #         if abs(reaction['flux']) <= ABS_TOL:
        #             reaction['direction'] = Direction.FORWARD
        #         else:
        #             assert reaction['flux'] < +ABS_TOL
        #             reaction['direction'] = Direction.REVERSIBLE
        #     elif reaction['direction'] == Direction.REVERSE:
        #         reaction['direction'] = Direction.REVERSIBLE
        #     else:
        #         assert False

    for i in range(len(data['reactions'])):
        forward = data['reactions'][i]
        if forward['direction'] is Direction.REVERSIBLE or forward['direction'] is Direction.DAMMED:
            continue

        forward_id = forward['id']
        if forward_id.endswith('_L2R'):
            reverse_id = '{}_R2L'.format(forward_id[: -4])
        elif forward_id.endswith('_R2L'):
            reverse_id = '{}_L2R'.format(forward_id[: -4])
        else:
            continue

        if reverse_id not in reaction_ids:
            continue

        reverse = data['reactions'][reaction_ids[reverse_id]]
        if reverse['direction'] is not Direction.DAMMED:
            log_.info('[{}, {}] is reversible now.'.format(forward_id, reverse_id))
            forward['direction'] = Direction.REVERSIBLE
            reverse['direction'] = Direction.REVERSIBLE

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
            metname = None
            for name in reaction['metabolites']:
                if re.sub(r'[-\+]', '_', rid).endswith('_{}'.format(name[: -2])):
                    if metname is None or len(metname) < len(name):
                        metname = name
            if metname is None:
                log_.warn('reaction [{}] has reactants [{}].'.format(reaction['name'].split('/')[0].replace('-', '_'), [name for name, coef in reaction['metabolites'].items()]))
                continue
            rid = rid[: -len(metname)+1]

        loc = default_location
        mobj = re.search(r'\[([^\[\]]+)\]', rid)
        if mobj is not None:
            loc = mobj[1]
            rid = rid[: -len(loc)-2]
        assert '[' not in rid and ']' not in rid

        reaction['source'] = dict(id=rid, location=loc)
        log_.debug('reaction id of [{}] was assumed to be [{}@{}].'.format(reaction['id'], rid, loc))

    return data

def check_reversibility(data):
    reactions = data['reactions']

    def compare_metabolites(m1, m2, abs_tol=1e-12, rel_tol=1e-6):
        if set(m1) != set(m2):
            return 0
        elif all(abs(float(val) + float(m2[key])) < abs_tol + rel_tol * abs(float(val)) for key, val in m1.items()):
            return -1
        elif all(abs(float(val) - float(m2[key])) < abs_tol + rel_tol * abs(float(val)) for key, val in m1.items()):
            return +1
        return 0

    done = [False] * len(reactions)
    fluxes = []
    for i, reaction1 in enumerate(reactions):
        if done[i]:
            continue
        fluxes.append([(i, +1)])
        for j_, reaction2 in enumerate(reactions[i+1: ]):
            j = i + 1 + j_
            if done[j]:
                continue
            res = compare_metabolites(reaction1['metabolites'], reaction2['metabolites'])
            if res == 0:
                pass
            else:
                done[j] = True
                fluxes[-1].append((j, res))

    return fluxes

def generate_ecocyc_fba(ECOCYC_VERSION="21.1", showall=False):
    log_.info('generate_ecocyc_fba(ECOCYC_VERSION="{}")'.format(ECOCYC_VERSION))

    from . import ecocyc
    ecocyc.load(path=INPUTS_PATH, version=ECOCYC_VERSION)

    filename = os.path.join(INPUTS_PATH, ECOCYC_VERSION, 'data/fba/fba-examples/ecocyc-21.0-gem-cs-glucose-tea-none.json')
    # filename = os.path.join(INPUTS_PATH, ECOCYC_VERSION, 'data/fba/fba-examples/ecocyc-21.0-gem-full-nutrient-set-cs-glucose-tea-oxygen.json')
    if not os.path.isfile(filename):
        raise RuntimeError('an input file [{}] could not be found'.format(filename))
    log_.info('read a file [{}]'.format(filename))
    with open(filename, 'r') as fin:
        data = json.load(fin)

    log_.info('[{}] reactions are found.'.format(len(data['reactions'])))
    # log_.info('[{}] reversible reactions are found.'.format(len([reaction for reaction in data['reactions'] if 'lower_bound' in reaction and 'upper_bound' in reaction and (reaction['lower_bound'] < 0 and reaction['upper_bound'] > 0)])))

    data = assume_source(data)
    data = solve_fba_using_cobra(data)

    for reaction in data['reactions']:
        if 'source' in reaction:
            r = ecocyc.find_reaction(reaction['source']['id'])
            if r is not None:
                log_.debug('{} <= {}'.format(reaction['id'], str(r)))
                if 'ENZYMATIC-REACTION' in r:
                    for enzrxn_id in r['ENZYMATIC-REACTION']:
                        enzrxn = ecocyc.find_enzrxn(enzrxn_id)
                        if 'ENZYME' in enzrxn:
                            if 'enzyme' not in reaction:
                                reaction['enzyme'] = [enzrxn['ENZYME']]
                            else:
                                reaction['enzyme'].append(enzrxn['ENZYME'])
                        # reaction['cofactor'].append(enzrxn['COFACTOR'])
            else:
                log_.warn('Reaction [{}] has no corresponding source named [{}] [{}]'.format(reaction['id'], reaction['source']['id'], reaction['name']))
        else:
            log_.debug('Reaction [{}] has no source'.format(reaction['id']))

    # fluxes = check_reversibility(data)
    # reactions = data['reactions']
    # for flux in fluxes:
    #     if len(flux) < 2:
    #         continue
    #     elif all(reactions[idx]['direction'] in (Direction.REVERSIBLE, Direction.DAMMED) for idx, coef in flux):
    #         continue
    #     # elif all(abs(reactions[idx]['flux']) <= ABS_TOL for idx, coef in flux):
    #     #     continue

    #     if all(reactions[idx]['direction'] in (Direction.REVERSIBLE, Direction.DAMMED) for idx, coef in flux):
    #         for idx, coef in flux:
    #             if reactions[idx]['direction'] not in (Direction.REVERSIBLE, Direction.DAMMED):
    #                 name = reactionss[idx]['id']
    #                 log_.warn('[{}] could be reversible.'.format(name))
    #                 reactions[idx]['direction'] = Direction.REVERSIBLE
    #     elif (any((reactions[idx]['direction'] is Direction.FORWARD and coef > 0)
    #               or (reactions[idx]['direction'] is Direction.REVERSE and coef < 0) 
    #               for idx, coef in flux)
    #           and any((reactions[idx]['direction'] is Direction.FORWARD and coef < 0)
    #               or (reactions[idx]['direction'] is Direction.REVERSE and coef > 0) 
    #               for idx, coef in flux)):
    #         for idx, coef in flux:
    #             name = reactions[idx]['id']
    #             log_.warn('[{}] could be reversible.'.format(name))
    #             reactions[idx]['direction'] = Direction.REVERSIBLE

    #     # tot = sum(reactions[idx]['flux'] * coef for idx, coef in flux)
    #     # forward = [(idx, coef) for idx, coef in flux
    #     #            if (reactions[idx]['direction'] is Direction.FORWARD and coef > 0) or (reactions[idx]['direction'] is Direction.REVERSE and coef < 0)]
    #     # reverse = [(idx, coef) for idx, coef in flux
    #     #            if (reactions[idx]['direction'] is Direction.FORWARD and coef < 0) or (reactions[idx]['direction'] is Direction.REVERSE and coef > 0)]
    #     # reversible = [(idx, coef) for idx, coef in flux if reactions[idx]['direction'] == Direction.REVERSIBLE]

    #     # if abs(tot) <= ABS_TOL:
    #     #     V = 1.0  # default velocity
    #     #     if len(forward) != 0 and len(reverse) != 0:
    #     #         forward.extend((idx, coef) for idx, coef in reversible if coef > 0)
    #     #         reverse.extend((idx, coef) for idx, coef in reversible if coef < 0)
    #     #     elif len(reversible) != 0 and (len(forward) != 0 or len(reverse) != 0):
    #     #         if len(forward) == 0:
    #     #             forward = reversible
    #     #         else:
    #     #             reverse = reversible
    #     #     else:
    #     #         continue  # do nothing
    #     #     for idx, coef in forward:
    #     #         reactions[idx]['flux'] = V / len(forward) * (+1 if coef > 0 else -1)
    #     #     for idx, coef in reverse:
    #     #         reactions[idx]['flux'] = V / len(reverse) * (+1 if coef < 0 else -1)
    #     # elif tot > 0.0:
    #     #     if len(forward) == 0:
    #     #         assert len(reversible) != 0
    #     #         forward = reversible
    #     #     else:
    #     #         forward.extend((idx, coef) for idx, coef in reversible if coef > 0)
    #     #         reverse.extend((idx, coef) for idx, coef in reversible if coef < 0)
    #     #     if len(reverse) == 0:
    #     #         for idx, coef in forward:
    #     #             reactions[idx]['flux'] = tot / len(forward) * (+1 if coef > 0 else -1)
    #     #     else:
    #     #         for idx, coef in forward:
    #     #             reactions[idx]['flux'] = 2 * tot / len(forward) * (+1 if coef > 0 else -1)
    #     #         for idx, coef in reverse:
    #     #             reactions[idx]['flux'] = tot / len(reverse) * (+1 if coef < 0 else -1)
    #     # else:  # tot < 0.0
    #     #     if len(reverse) == 0:
    #     #         assert len(reversible) != 0
    #     #         reverse = reversible
    #     #     else:
    #     #         forward.extend((idx, coef) for idx, coef in reversible if coef > 0)
    #     #         reverse.extend((idx, coef) for idx, coef in reversible if coef < 0)
    #     #     if len(forward) == 0:
    #     #         for idx, coef in reverse:
    #     #             reactions[idx]['flux'] = -tot / len(forward) * (+1 if coef < 0 else -1)
    #     #     else:
    #     #         for idx, coef in forward:
    #     #             reactions[idx]['flux'] = -tot / len(forward) * (+1 if coef > 0 else -1)
    #     #         for idx, coef in reverse:
    #     #             reactions[idx]['flux'] = -2 * tot / len(reverse) * (+1 if coef < 0 else -1)

    compounds = []
    for reaction in data['reactions']:
        flux = reaction['flux']
        if not showall and (abs(flux) <= ABS_TOL and reaction['direction'] is not Direction.REVERSIBLE):
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

    def calc_params(flux, metabolites, reversible=Direction.REVERSIBLE):
        if abs(flux) <= ABS_TOL: # and reversible is not Direction.REVERSIBLE:
            return (0.0, 0.0)

        Km = 1.0
        denom1 = functools.reduce(lambda x, y: (x * pow(1 + compounds[y[0]] / Km, int(-y[1]))) if y[1] < 0 else x, metabolites.items(), 1.0)
        denom1 += functools.reduce(lambda x, y: (x * pow(1 + compounds[y[0]] / Km, int(y[1]))) if y[1] > 0 else x, metabolites.items(), 1.0)
        denom1 -= 1

        if reversible is Direction.REVERSIBLE or reversible is Direction.DAMMED:
            Keq = 2.0 if flux > 0 else 0.5
            num = functools.reduce(lambda x, y: (x * pow(compounds[y[0]], int(-y[1]))) if y[1] < 0 else x, metabolites.items(), 1.0)
            num -= functools.reduce(lambda x, y: (x * pow(compounds[y[0]], int(y[1]))) if y[1] > 0 else x, metabolites.items(), 1.0) / Keq
            if abs(flux) <= ABS_TOL:
                Vmax = 1.0 / (num / denom1)
                return (Vmax, Vmax)
            else:
                Vmax = flux / (num / denom1)
                return (Vmax, Vmax / Keq)
        elif reversible is Direction.FORWARD:
            assert flux > ABS_TOL
            num = functools.reduce(lambda x, y: (x * pow(compounds[y[0]], int(-y[1]))) if y[1] < 0 else x, metabolites.items(), 1.0)
            Vmax = flux / (num / denom1)
            return (Vmax, 0.0)
        else:
            assert reversible is Direction.REVERSE
            assert flux < -ABS_TOL
            num = -functools.reduce(lambda x, y: (x * pow(compounds[y[0]], int(y[1]))) if y[1] > 0 else x, metabolites.items(), 1.0)
            Vmax = flux / (num / denom1)
            return (0.0, Vmax)

    filename = os.path.join(OUTPUTS_PATH, 'metabolism.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')
        for reaction in data['reactions']:
            flux = reaction['flux']
            if not showall and abs(flux) <= ABS_TOL:
            # if not showall and (abs(flux) <= ABS_TOL and reaction['direction'] is not Direction.REVERSIBLE):
                continue
            assert ';' not in reaction['id'] and ':' not in reaction['id']
            assert all(';' not in name and ':' not in name for name in reaction['metabolites'])

            # velocity = lambda x: x / (x + 0.25)
            velocity = lambda x: x

            if reaction['id'] == 'MetaFlux_obj':
                assert abs(flux) > ABS_TOL
                for name, coef in reaction['metabolites'].items():
                    if coef < 0:
                        vfmax, vrmax = calc_params(flux * -coef, {name: -1.0}, Direction.FORWARD)
                        writer.writerow(('{}_{}'.format(reaction['id'], name), '{}:{}'.format(name, 1.0), '', vfmax, vrmax, ''))
                    elif coef > 0:
                        vfmax, vrmax = calc_params(flux * coef, {name: +1.0}, Direction.FORWARD)
                        writer.writerow(('{}_{}'.format(reaction['id'], name), '', '{}:{}'.format(name, 1.0), vfmax, vrmax, ''))
                    else:
                        assert False  # never get here
                continue

            for name, coef in reaction['metabolites'].items():
                coef = abs(coef)
                if coef - floor(coef) > ABS_TOL:
                    print(reaction)
                    raise RuntimeError('[{}] has non-int coef [{}].'.format(name, coef))

            reactants = ';'.join('{}:{}'.format(name, -coef) for name, coef in reaction['metabolites'].items() if coef < 0)
            products = ';'.join('{}:{}'.format(name, coef) for name, coef in reaction['metabolites'].items() if coef > 0)

            # vfmax, vrmax = calc_params(flux, reaction['metabolites'], reaction['direction'])
            vfmax, vrmax = calc_params(flux, reaction['metabolites'], Direction.REVERSIBLE)

            writer.writerow((reaction['id'], reactants, products, vfmax, vrmax, ';'.join(reaction.get('enzyme', ()))))

def generate_ecocyc_plexes(ECOCYC_VERSION="21.1"):
    log_.info('generate_ecocyc_fba(ECOCYC_VERSION="{}")'.format(ECOCYC_VERSION))

    from . import ecocyc
    ecocyc.load(path=INPUTS_PATH, version=ECOCYC_VERSION)

    filename = os.path.join(OUTPUTS_PATH, 'plexes.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')
        for feature in ecocyc.proteins():
            if 'UNMODIFIED-FORM' in feature:
                writer.writerow((feature['UNIQUE-ID'], '{}:1'.format(feature['UNMODIFIED-FORM'])))
            elif 'COMPONENTS' in feature:
                writer.writerow((feature['UNIQUE-ID'], ';'.join('{}:{:g}'.format(component, coef) for component, coef in feature['COMPONENTS'])))
            else:
                continue

def generate_ecocyc_gene_product_map(ECOCYC_VERSION="21.1"):
    log_.info('generate_ecocyc_gene_product_map(ECOCYC_VERSION="{}")'.format(ECOCYC_VERSION))

    from . import ecocyc
    ecocyc.load(path=INPUTS_PATH, version=ECOCYC_VERSION)

    filename = os.path.join(OUTPUTS_PATH, 'gene_product_map.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')

        for feature in ecocyc.proteins():
            if 'UNMODIFIED-FORM' not in feature and 'COMPONENTS' not in feature and 'GENE' in feature:
                writer.writerow((feature['GENE'], feature['UNIQUE-ID']))
        for feature in ecocyc.rnas():
            if 'UNMODIFIED-FORM' not in feature and 'GENE' in feature:
                writer.writerow((feature['GENE'], feature['UNIQUE-ID']))

def generate_ecocyc_transunits(ECOCYC_VERSION="21.1"):
    log_.info('generate_ecocyc_transunits(ECOCYC_VERSION="{}")'.format(ECOCYC_VERSION))

    from . import ecocyc
    ecocyc.load(path=INPUTS_PATH, version=ECOCYC_VERSION)

    transunits = defaultdict(list)
    for feature in ecocyc.transunits():
        if 'COMPONENTS' in feature:
            for component in feature['COMPONENTS']:
                transunits[component].append(feature['UNIQUE-ID'])

    filename = os.path.join(OUTPUTS_PATH, 'transunits.csv')
    with open(filename, 'w') as fout:
        log_.info('output a file [{}]'.format(filename))
        writer = csv.writer(fout, lineterminator='\n')

        for feature in ecocyc.proteins():
            if 'UNMODIFIED-FORM' not in feature and 'COMPONENTS' not in feature and 'GENE' in feature:
                for transunit in transunits[feature['GENE']]:
                    writer.writerow((transunit, feature['UNIQUE-ID']))
        # for feature in ecocyc.rnas():
        #     if 'UNMODIFIED-FORM' not in feature and 'GENE' in feature:
        #         writer.writerow((feature['GENE'], feature['UNIQUE-ID']))


if __name__ == "__main__":
    basicConfig(level=INFO)

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

    generate_ecocyc_fba(showall=False)
    generate_ecocyc_plexes()
    generate_ecocyc_gene_product_map()
    generate_ecocyc_transunits()