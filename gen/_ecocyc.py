#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import copy
import os.path
import itertools
import logging

from . import utils

__program__ = '_ecocyc.py'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu'

ECOCYC_VERSION = '21.1'
ECOCYC_DATA_PATH = os.path.join('.', '{}/data'.format(ECOCYC_VERSION))
# ECOCYC_DATA_PATH = os.path.join(utils.DATA_PATH, '{}/data'.format(ECOCYC_VERSION))
ECOCYC_PROTEINS = []
ECOCYC_COMPOUNDS = []
ECOCYC_CLASSES = []
ECOCYC_RNAS = []
ECOCYC_GENES = []
ECOCYC_TRANSUNITS = []
ECOCYC_TERMINATORS = []
ECOCYC_PROMOTERS = []
ECOCYC_DNABINDSITES = []
ECOCYC_REACTIONS = []
ECOCYC_ENZRXNS = []
ECOCYC_BINDRXNS = []
ECOCYC_REGULONS = []
ECOCYC_PROTLIGANDCPLXES = []
ECOCYC_PATHWAYS = []
ECOCYC_PROTEIN_FEATURES = []
# ECOCYC_CMR_GENE_LOOKUP = []
# ECOCYC_GENE_LINKS = []

__all__ = ['load', 'ECOCYC_VERSION',
           'ECOCYC_PROTEINS', 'ECOCYC_COMPOUNDS', 'ECOCYC_CLASSES',
           'ECOCYC_RNAS', 'ECOCYC_GENES', 'ECOCYC_TRANSUNITS',
           'ECOCYC_TERMINATORS', 'ECOCYC_PROMOTERS', 'ECOCYC_DNABINDSITES',
           'ECOCYC_REACTIONS', 'ECOCYC_ENZRXNS', 'ECOCYC_BINDRXNS',
           'ECOCYC_REGULONS', 'ECOCYC_PROTLIGANDCPLXES', 'ECOCYC_PATHWAYS',
           'ECOCYC_PROTEIN_FEATURES']

logger = logging.getLogger('libecoli')


def read_ecocyc_file(filename, encoding='ISO-8859-1'):
    def split_line(line):
        data = line.split(' - ')
        if len(data) > 1:
            return (data[0], ' - '.join(data[1:]))
        elif line[-2:] == ' -':
            return (line[: -2], '')
        else:
            raise RuntimeError("invalid format: " + repr(data))

    attributes = []
    unique_id, cache = None, []
    retval = {}
    linecnt = 1
    with open(filename, 'r', encoding=encoding) as fin:
        while True:
            line = fin.readline()
            linecnt += 1
            if line is None or line == '':
                break
            line = line.strip()
            if line[0] == '#':
                pass
            elif line == '//':
                if unique_id is None:
                    raise RuntimeError('No UNIQUE-ID was specified.')
                retval[unique_id] = tuple(cache)
                unique_id, cache = None, []
            elif line[0] == '/':
                cache[-1] = cache[-1][: -1] + (cache[-1][-1] + line[1:],)
            elif line[0] == '^':
                key, value = split_line(line[1:])
                attributes.append(key)
                cache[-1] = cache[-1] + (key, value)
            else:
                try:
                    key, value = split_line(line)
                except RuntimeError as err:
                    # XXX: exception for 'pathways.dat'
                    if len(cache[-1]) == 2 and cache[-1][0] == 'COMMENT':
                        logger.warning(
                            "invalid format in 'COMMENT' of '{0}'".format(
                                unique_id)
                            + ". '/' seems missing before '{0}': ".format(line)
                            + filename)
                        cache[-1] = cache[-1][: -1] + (cache[-1][-1] + line,)
                        continue
                    raise RuntimeError('At Line {}: {}'.format(linecnt, str(err)))

                attributes.append(key)
                if key == 'UNIQUE-ID':
                    unique_id = value
                cache.append((key, value))
    logger.debug(
        "{0} {1:s} {2:d}".format(filename, str(sorted(set(attributes))), len(retval)))
    return retval


def get_attributes(entry, key):
    return filter(lambda attr: attr[0] == key, entry)


def has_key(entry, key):
    return any([attr[0] == key for attr in entry])


def get_value(entry, key):
    retval = []
    for attr in get_attributes(entry, key):
        if len(attr) != 2:
            raise ValueError(
                'length of a value must be two.'
                + ' {0:d} was given.'.format(len(attr)))
        retval.append(attr[1])
    return tuple(retval)


def distinct(entries, key):
    return set(itertools.chain(
        *[[attr[1] for attr in entry if attr[0] == key] for entry in entries]))


def set_value(entry, key, target, unique=False):
    if not has_key(entry, key):
        return False
    if not unique:
        target[key] = get_value(entry, key)
    else:
        v = get_value(entry, key)
        if len(v) != 1:
            raise ValueError('value must have only one element: {}'.format(repr(v)))
        target[key] = v[0]
    return True


def unwrap(name):
    if len(name) > 1 and name[0] == '|' and name[-1] == '|':
        return name[1: -1]
    return name


def load_proteins():
    """'ABBREV-NAME', 'ATOM-CHARGES', 'CATALYZES', 'CHEMICAL-FORMULA',
    'CITATIONS', 'COEFFICIENT', 'COFACTORS-OF', 'COMMENT', 'COMMENT-INTERNAL',
    'COMMON-NAME', 'COMPONENT-COEFFICIENTS', 'COMPONENT-OF', 'COMPONENTS',
    'CONSENSUS-SEQUENCE', 'CREDITS', 'DBLINKS', 'DNA-FOOTPRINT-SIZE',
    'ENZYME-NOT-USED-IN', 'FEATURES', 'GENE', 'GO-TERMS', 'IN-MIXTURE',
    'INCHI', 'ISOZYME-SEQUENCE-SIMILARITY', 'LOCATIONS', 'MODIFIED-FORM',
    'MOLECULAR-WEIGHT', 'MOLECULAR-WEIGHT-EXP', 'MOLECULAR-WEIGHT-KD',
    'MOLECULAR-WEIGHT-SEQ', 'MONOISOTOPIC-MW', 'NEIDHARDT-SPOT-NUMBER',
    'NON-STANDARD-INCHI', 'PI', 'PROMOTER-BOX-NAME-1', 'PROMOTER-BOX-NAME-2',
    'RECOGNIZED-PROMOTERS', 'REGULATED-BY', 'REGULATES', 'SMILES', 'SPECIES',
    'SPLICE-FORM-INTRONS', 'STRUCTURE-BONDS', 'SYMMETRY', 'SYNONYMS',
    'TEMPLATE-FILE', 'TYPES', 'UNIQUE-ID', 'UNMODIFIED-FORM'
    """
    global ECOCYC_PROTEINS

    filename = os.path.join(ECOCYC_DATA_PATH, 'proteins.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_PROTEINS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'GENE', value, True)
        set_value(entry, 'TYPES', value)
        set_value(entry, 'COMMON-NAME', value, True)
        set_value(entry, 'ABBREV-NAME', value, True)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'FEATURES', value)
        set_value(entry, 'REGULATES', value)

        set_value(entry, 'INCHI', value, True)
        set_value(entry, 'NON-STANDARD-INCHI', value, True)

        if set_value(entry,'DNA-FOOTPRINT-SIZE', value):
            value['DNA-FOOTPRINT-SIZE'] = [
                int(footprint) for footprint in value['DNA-FOOTPRINT-SIZE']]

        if has_key(entry, 'COMPONENTS'):
            value['COMPONENTS'] = []
            for attr in get_attributes(entry, 'COMPONENTS'):
                if len(attr) == 2:
                    value['COMPONENTS'].append((attr[1], 1.0))
                elif len(attr) == 4 and attr[2] == 'COEFFICIENT':
                    value['COMPONENTS'].append((attr[1], float(attr[3])))
                else:
                    raise ValueError(
                        'invalid format in COMPONENTS: ' + repr(value))
            value['COMPONENTS'] = tuple(value['COMPONENTS'])

            for target, coef in value['COMPONENTS']:
                if target in data.keys():
                    if (not has_key(data[target], 'COMPONENT-OF')
                            or unique_id not in get_value(
                                data[target], 'COMPONENT-OF')):
                        raise ValueError(
                            repr(entry) + ' not in ' + repr(data[target]))

        if set_value(entry, 'UNMODIFIED-FORM', value, True):
            target = value['UNMODIFIED-FORM']
            if target in data.keys():
                if (not has_key(data[target], 'MODIFIED-FORM')
                        or unique_id not in get_value(
                            data[target], 'MODIFIED-FORM')):
                    raise ValueError(
                        repr(entry) + ' not in ' + repr(data[target]))

        ECOCYC_PROTEINS.append(value)


def load_compounds():
    """'ABBREV-NAME', 'ATOM-CHARGES', 'CFG-ICON-COLOR', 'CHEMICAL-FORMULA',
    'CITATIONS', 'COFACTORS-OF', 'COMMENT', 'COMMENT-INTERNAL', 'COMMON-NAME',
    'COMPONENT-OF', 'COMPONENTS', 'CREDITS', 'DBLINKS', 'GROUP-COORDS-2D',
    'GROUP-INTERNALS', 'HAS-NO-STRUCTURE?', 'IN-MIXTURE', 'INCHI',
    'INTERNALS-OF-GROUP', 'MOLECULAR-WEIGHT', 'MONOISOTOPIC-MW', 'N+1-NAME',
    'N-NAME', 'NON-STANDARD-INCHI', 'PKA1', 'PKA2', 'PKA3', 'RADICAL-ATOMS',
    'REGULATES', 'SMILES', 'STRUCTURE-GROUPS', 'STRUCTURE-LINKS', 'SUPERATOMS',
    'SYNONYMS', 'SYSTEMATIC-NAME', 'TAUTOMERS', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_COMPOUNDS

    filename = os.path.join(ECOCYC_DATA_PATH, 'compounds.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_COMPOUNDS = []
    for unique_id, entry in data.items():
        value = {}
        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'COMMON-NAME', value, True)
        set_value(entry, 'ABBREV-NAME', value, True)
        set_value(entry, 'SYNONYMS', value)
        # set_value(entry, 'REGULATES', value)

        set_value(entry, 'DBLINKS', value)
        set_value(entry, 'INCHI', value, True)
        set_value(entry, 'NON-STANDARD-INCHI', value, True)

        set_value(entry, 'SMILES', value, True)
        set_value(entry, 'COMPONENTS', value)
        if set_value(entry, 'CHEMICAL-FORMULA', value):
            formula = {}
            for target in value['CHEMICAL-FORMULA']:
                atom, count = target[1: -1].split(' ')
                formula[atom] = int(count)
            value['CHEMICAL-FORMULA'] = formula

        ECOCYC_COMPOUNDS.append(value)


def load_rnas():
    """'ANTICODON', 'CITATIONS', 'COMMENT', 'COMMENT-INTERNAL', 'COMMON-NAME',
    'COMPONENT-OF', 'CREDITS', 'DBLINKS', 'GENE', 'GO-TERMS', 'LOCATIONS',
    'MODIFIED-FORM', 'REGULATES', 'SPECIES', 'SYNONYMS', 'TYPES', 'UNIQUE-ID',
    'UNMODIFIED-FORM'
    """
    global ECOCYC_RNAS

    filename = os.path.join(ECOCYC_DATA_PATH, 'rnas.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_RNAS = []
    for unique_id, entry in data.items():
        value = {}
        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'GENE', value, True)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'COMMON-NAME', value, True)
        # set_value(entry, 'REGULATES', value)

        # if unique_id == 'charged-ileX-tRNA':  # XXX: EXCEPTION
        #     v = get_value(entry, 'UNMODIFIED-FORM')
        #     logging.warning(
        #         unique_id + ' has multiple UNMODIFIED-FORM: ' + repr(v))
        #     value['UNMODIFIED-FORM'] = v[-1]
        # else:
        #     set_value(entry, 'UNMODIFIED-FORM', value, True)
        set_value(entry, 'UNMODIFIED-FORM', value)

        ECOCYC_RNAS.append(value)


def load_classes():
    """'COMMENT', 'COMMON-NAME', 'SYNONYMS', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_CLASSES

    filename = os.path.join(ECOCYC_DATA_PATH, 'classes.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_CLASSES = []
    for unique_id, entry in data.items():
        value = {}
        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'COMMON-NAME', value, True)
        set_value(entry, 'TYPES', value)
        ECOCYC_CLASSES.append(value)


def load_genes():
    """'ACCESSION-1', 'ACCESSION-2', 'CENTISOME-POSITION', 'CITATIONS',
    'COMMENT', 'COMMENT-INTERNAL', 'COMMON-NAME', 'COMPONENT-OF', 'DBLINKS',
    'IN-GROUP', 'INTERRUPTED?', 'KNOCKOUT-GROWTH-OBSERVATIONS', 'LAST-UPDATE',
    'LEFT-END-POSITION', 'MEMBER-SORT-FN', 'PRODUCT', 'REGULATED-BY',
    'RIGHT-END-POSITION', 'SYNONYMS', 'TRANSCRIPTION-DIRECTION', 'TYPES',
    'UNIQUE-ID'
    """
    global ECOCYC_GENES

    filename = os.path.join(ECOCYC_DATA_PATH, 'genes.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_GENES = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'COMMON-NAME', value, True)
        set_value(entry, 'TYPES', value)

        # assert(set_value(entry, 'COMPONENT-OF', value))
        if not set_value(entry, 'COMPONENT-OF', value):
            #XXX: rather use 'UNMAPPED-COMPONENT-OF'
            assert value['UNIQUE-ID'] in ('G0-9281', 'G26')
            logging.warning(
                unique_id + " has no 'COMPONENT-OF': " + repr(entry))
            value['COMPONENT-OF'] = ('COLI-K12', )
        set_value(entry, 'PRODUCT', value)
        # set_value(entry, 'REGULATED-BY', value)

        set_value(entry, 'ACCESSION-1', value, True)  # b-number?
        # set_value(entry, 'ACCESSION-2', value, True)  # EcoliWiki

        if set_value(entry, 'LEFT-END-POSITION', value, True):
            value['LEFT-END-POSITION'] = int(value['LEFT-END-POSITION'])
        if set_value(entry, 'RIGHT-END-POSITION', value, True):
            value['RIGHT-END-POSITION'] = int(value['RIGHT-END-POSITION'])

        if set_value(entry, 'TRANSCRIPTION-DIRECTION', value, True):
            if value['TRANSCRIPTION-DIRECTION'] not in ('+', '-'):
                raise ValueError(
                    'invalid direction was given: '
                    + value['TRANSCRIPTION-DIRECTION'])
            value['TRANSCRIPTION-DIRECTION'] = (
                value['TRANSCRIPTION-DIRECTION'] == '+')

        ECOCYC_GENES.append(value)


def validate_genes():
    global ECOCYC_GENES

    parents = [entry['UNIQUE-ID'] for entry in ECOCYC_TRANSUNITS]

    for entry in ECOCYC_GENES:
        assert('UNIQUE-ID' in entry.keys())
        assert('COMPONENT-OF' in entry.keys())

        for target in entry['COMPONENT-OF']:
            assert(isinstance(target, str))
            assert(target in parents or target.startswith('COLI-K12'))

            if not target.startswith('COLI-K12'):
                parent = ECOCYC_TRANSUNITS[parents.index(target)]
                assert('COMPONENTS' in parent.keys())
                assert(entry['UNIQUE-ID'] in parent['COMPONENTS'])

        if 'ACCESSION-1' in entry.keys():
            assert(entry['ACCESSION-1'][0] == 'b')

        if 'TRANSCRIPTION-DIRECTION' in entry.keys():
            assert(isinstance(entry['TRANSCRIPTION-DIRECTION'], bool))


def load_transunits():
    """'CITATIONS', 'COMMENT', 'COMMON-NAME', 'COMPONENTS', 'EXTENT-UNKNOWN?',
    'REGULATED-BY', 'SYNONYMS', 'TRANSCRIPTION-DIRECTION', 'TYPES',
    'UNIQUE-ID'
    """
    global ECOCYC_TRANSUNITS

    filename = os.path.join(ECOCYC_DATA_PATH, 'transunits.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_TRANSUNITS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'COMMON-NAME', value, True)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'TYPES', value)

        assert(set_value(entry, 'COMPONENTS', value))
        # set_value(entry, 'REGULATED-BY', value)

        if set_value(entry, 'TRANSCRIPTION-DIRECTION', value, True):
            if value['TRANSCRIPTION-DIRECTION'] not in ('+', '-'):
                raise ValueError(
                    'invalid direction was given: '
                    + value['TRANSCRIPTION-DIRECTION'])
            value['TRANSCRIPTION-DIRECTION'] = (
                value['TRANSCRIPTION-DIRECTION'] == '+')

        ECOCYC_TRANSUNITS.append(value)


def validate_transunits():
    global ECOCYC_TRANSUNITS

    # logger.debug(ECOCYC_TRANSUNITS)

    elements = [entry['UNIQUE-ID'] for entry in itertools.chain(
        ECOCYC_PROMOTERS, ECOCYC_TERMINATORS, ECOCYC_GENES,
        ECOCYC_DNABINDSITES)]

    for entry in ECOCYC_TRANSUNITS:
        assert('UNIQUE-ID' in entry.keys())
        assert('COMPONENTS' in entry.keys())

        for target in entry['COMPONENTS']:
            assert(isinstance(target, str))

            if (target not in elements):
                logger.warning(
                    "'{0}' in 'COMPONENTS' of '{1}' is undefined".format(
                        target, entry['UNIQUE-ID']))

        if 'TRANSCRIPTION-DIRECTION' in entry.keys():
            assert(isinstance(entry['TRANSCRIPTION-DIRECTION'], bool))


def load_terminators():
    """'CITATIONS', 'COMMENT', 'COMPONENT-OF', 'LEFT-END-POSITION',
    'REGULATED-BY', 'RIGHT-END-POSITION', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_TERMINATORS

    filename = os.path.join(ECOCYC_DATA_PATH, 'terminators.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_TERMINATORS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)

        if set_value(entry, 'LEFT-END-POSITION', value, True):
            value['LEFT-END-POSITION'] = int(value['LEFT-END-POSITION'])
        if set_value(entry, 'RIGHT-END-POSITION', value, True):
            value['RIGHT-END-POSITION'] = int(value['RIGHT-END-POSITION'])

        ECOCYC_TERMINATORS.append(value)


def load_promoters():
    """'ABSOLUTE-PLUS-1-POS', 'BINDS-SIGMA-FACTOR', 'CITATIONS', 'COMMENT',
    'COMMENT-INTERNAL', 'COMMON-NAME', 'COMPONENT-OF', 'MINUS-10-LEFT',
    'MINUS-10-RIGHT', 'MINUS-35-LEFT', 'MINUS-35-RIGHT', 'REGULATED-BY',
    'SYNONYMS', 'TRANSCRIPTION-DIRECTION', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_PROMOTERS

    filename = os.path.join(ECOCYC_DATA_PATH, 'promoters.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_PROMOTERS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'SYNONYMS', value)

        if set_value(entry, 'ABSOLUTE-PLUS-1-POS', value, True):
            value['ABSOLUTE-PLUS-1-POS'] = int(value['ABSOLUTE-PLUS-1-POS'])
        if set_value(entry, 'MINUS-10-LEFT', value, True):
            value['MINUS-10-LEFT'] = int(value['MINUS-10-LEFT'])
        if set_value(entry, 'MINUS-10-RIGHT', value, True):
            value['MINUS-10-RIGHT'] = int(value['MINUS-10-RIGHT'])
        if set_value(entry, 'MINUS-35-LEFT', value, True):
            value['MINUS-35-LEFT'] = int(value['MINUS-35-LEFT'])
        if set_value(entry, 'MINUS-35-RIGHT', value, True):
            value['MINUS-35-RIGHT'] = int(value['MINUS-35-RIGHT'])

        # set_value(entry, 'REGULATED-BY', value)
        # set_value(entry, 'COMPONENT-OF', value)

        set_value(entry, 'BINDS-SIGMA-FACTOR', value)  #XXX: This is deprecated in the latest version

        if set_value(entry, 'TRANSCRIPTION-DIRECTION', value, True):
            if value['TRANSCRIPTION-DIRECTION'] not in ('+', '-'):
                raise ValueError(
                    'invalid direction was given: '
                    + value['TRANSCRIPTION-DIRECTION'])
            value['TRANSCRIPTION-DIRECTION'] = (
                value['TRANSCRIPTION-DIRECTION'] == '+')

        ECOCYC_PROMOTERS.append(value)


def load_dnabindsites():
    """'ABS-CENTER-POS', 'CITATIONS', 'COMMENT', 'COMPONENT-OF', 'DBLINKS',
    'INVOLVED-IN-REGULATION', 'SITE-LENGTH', 'SYNONYMS', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_DNABINDSITES

    filename = os.path.join(ECOCYC_DATA_PATH, 'dnabindsites.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_DNABINDSITES = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'SYNONYMS', value)

        # set_value(entry, 'COMPONENT-OF', value)
        # set_value(entry, 'INVOLVED-IN-REGULATION', value)

        if set_value(entry, 'ABS-CENTER-POS', value, True):
            if value['ABS-CENTER-POS'] != 'NIL':
                value['ABS-CENTER-POS'] = float(value['ABS-CENTER-POS'])
            else:
                logger.warning(
                    '"NIL" was given for ABS-CENTER-POS of "'
                    + value['UNIQUE-ID']
                    + '" in dnabindsites.dat. Just ignored.')
                del value['ABS-CENTER-POS']
        if set_value(entry, 'SITE-LENGTH', value, True):
            value['SITE-LENGTH'] = int(value['SITE-LENGTH'])

        ECOCYC_DNABINDSITES.append(value)


def load_reactions():
    """'ATOM-MAPPINGS', 'CANNOT-BALANCE?', 'CITATIONS', 'COEFFICIENT',
    'COMMENT', 'COMMENT-INTERNAL', 'COMMON-NAME', 'COMPARTMENT',
    'CREDITS', 'DELTAG0', 'EC-NUMBER', 'ENZYMATIC-REACTION',
    'ENZYMES-NOT-USED', 'IN-PATHWAY', 'LEFT', 'MEMBER-SORT-FN', 'OFFICIAL?',
    'ORPHAN?', 'PHYSIOLOGICALLY-RELEVANT?', 'PREDECESSORS', 'PRIMARIES',
    'REACTION-DIRECTION', 'REACTION-LIST', 'REGULATED-BY', 'REQUIREMENTS',
    'RIGHT', 'RXN-LOCATIONS', 'SIGNAL', 'SPECIES', 'SPONTANEOUS?',
    'STD-REDUCTION-POTENTIAL', 'SYNONYMS', 'SYSTEMATIC-NAME',
    'TEMPLATE-FILE', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_REACTIONS

    filename = os.path.join(ECOCYC_DATA_PATH, 'reactions.dat')
    data = read_ecocyc_file(filename)

    tmp = []

    ECOCYC_REACTIONS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'COMMON-NAME', value)
        set_value(entry, 'SYNONYMS', value)

        if has_key(entry, 'EC-NUMBER'):
            # 'EC-NUMBER' may contain 'OFFICIAL?'
            value['EC-NUMBER'] = tuple(
                [v[1] for v in get_attributes(entry, 'EC-NUMBER')])

        # set_value(entry, 'REGULATED-BY', value)
        set_value(entry, 'ENZYMATIC-REACTION', value)
        set_value(entry, 'RXN-LOCATIONS', value)
        set_value(entry, 'IN-PATHWAY', value)

        if set_value(entry, 'REACTION-DIRECTION', value, True):
            assert(value['REACTION-DIRECTION'] in (
                'LEFT-TO-RIGHT', 'RIGHT-TO-LEFT', 'REVERSIBLE',
                'PHYSIOL-LEFT-TO-RIGHT', 'PHYSIOL-RIGHT-TO-LEFT',
                'IRREVERSIBLE-LEFT-TO-RIGHT'))

        for name in ('LEFT', 'RIGHT'):
            if has_key(entry, name):
                value[name] = []
                for attr in get_attributes(entry, name):
                    assert(len(attr) % 2 == 0)
                    v = {}
                    v['UNIQUE-ID'] = attr[1]
                    for i in range(1, len(attr) // 2):
                        attrname, attrvalue = attr[i * 2], attr[i * 2 + 1]
                        assert(attrname in ('COMPARTMENT', 'COEFFICIENT'))
                        v[attrname] = attrvalue
                    if 'COEFFICIENT' not in v.keys():
                        v['COEFFICIENT'] = '1'  #XXX: optional
                    value[name].append(v.copy())
                value[name] = tuple(value[name])

        ECOCYC_REACTIONS.append(value)


def load_enzrxns():
    """'ALTERNATIVE-COFACTORS', 'ALTERNATIVE-SUBSTRATES', 'CITATIONS',
    'COFACTOR-BINDING-COMMENT', 'COFACTORS', 'COMMENT', 'COMMENT-INTERNAL',
    'COMMON-NAME', 'CREDITS', 'ENZYME', 'KCAT', 'KM', 'PH-OPT',
    'PHYSIOLOGICALLY-RELEVANT?', 'REACTION', 'REACTION-DIRECTION',
    'REGULATED-BY', 'REQUIRED-PROTEIN-COMPLEX', 'SPECIFIC-ACTIVITY',
    'SYNONYMS', 'TEMPERATURE-OPT', 'TEMPLATE-FILE', 'TYPES', 'UNIQUE-ID',
    'VMAX'
    """
    global ECOCYC_ENZRXNS

    filename = os.path.join(ECOCYC_DATA_PATH, 'enzrxns.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_ENZRXNS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'COMMON-NAME', value)
        set_value(entry, 'SYNONYMS', value)

        set_value(entry, 'COFACTORS', value)
        set_value(entry, 'ENZYME', value, True)
        set_value(entry, 'REACTION', value)
        set_value(entry, 'REQUIRED-PROTEIN-COMPLEX', value)

        set_value(entry, 'SPECIFIC-ACTIVITY', value)

        if set_value(entry, 'REACTION-DIRECTION', value, True):
            assert(value['REACTION-DIRECTION'] in (
                'REVERSIBLE',
                'PHYSIOL-RIGHT-TO-LEFT', 'PHYSIOL-LEFT-TO-RIGHT',
                'IRREVERSIBLE-LEFT-TO-RIGHT', 'IRREVERSIBLE-RIGHT-TO-LEFT'))

        # set_value(entry, 'ALTERNATIVE-COFACTORS', value)  # XXX: format?
        # set_value(entry, 'ALTERNATIVE-SUBSTRATES', value)  # XXX: format?

        # set_value(entry, 'KCAT', value)
        # set_value(entry, 'KM', value)
        # set_value(entry, 'VMAX', value)
        # set_value(entry, 'REGULATED-BY', value)

        ECOCYC_ENZRXNS.append(value)


def load_bindrxns():
    """'ACTIVATORS', 'BASAL-TRANSCRIPTION-VALUE', 'CITATIONS', 'DBLINKS',
    'INHIBITORS', 'OFFICIAL-EC?', 'REACTANTS', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_BINDRXNS

    filename = os.path.join(ECOCYC_DATA_PATH, 'bindrxns.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_BINDRXNS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)

        for name in ('ACTIVATORS', 'INHIBITORS'):
            if set_value(entry, name, value):
                for v in value[name]:
                    assert(len(v.split('++')) == 2)
                value[name] = tuple(
                    [tuple(v.split('++')) for v in value[name]])

        set_value(entry, 'REACTANTS', value)

        ECOCYC_BINDRXNS.append(value)


def load_regulations():
    """'ACCESSORY-PROTEINS', 'ASSOCIATED-BINDING-SITE', 'CITATIONS',
    'COMMENT', 'CREDITS', 'DOWNSTREAM-GENES-ONLY?', 'GROWTH-CONDITIONS',
    'KI', 'MECHANISM', 'MODE', 'PAUSE-END-POS', 'PAUSE-START-POS',
    'PHYSIOLOGICALLY-RELEVANT?', 'REGULATED-ENTITY', 'REGULATOR',
    'SYNONYMS', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_REGULATIONS

    filename = os.path.join(ECOCYC_DATA_PATH, 'regulation.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_REGULATIONS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'MODE', value, True)

        set_value(entry, 'REGULATOR', value)
        set_value(entry, 'REGULATED-ENTITY', value, True)

        set_value(entry, 'ASSOCIATED-BINDING-SITE', value, True)
        set_value(entry, 'ACCESSORY-PROTEINS', value, True)

        ECOCYC_REGULATIONS.append(value)


def load_regulons():
    """'ABBREV-NAME', 'CATALYZES', 'CITATIONS', 'COMMENT', 'COMMENT-INTERNAL',
    'COMMON-NAME', 'COMPONENT-OF', 'COMPONENTS', 'CONSENSUS-SEQUENCE',
    'CREDITS', 'DBLINKS', 'DNA-FOOTPRINT-SIZE', 'FEATURES', 'GENE',
    'GO-TERMS', 'LOCATIONS', 'MODIFIED-FORM', 'MOLECULAR-WEIGHT-EXP',
    'MOLECULAR-WEIGHT-KD', 'MOLECULAR-WEIGHT-SEQ', 'PI', 'REGULATES',
    'SPECIES', 'SYMMETRY', 'SYNONYMS', 'TEMPLATE-FILE', 'TYPES',
    'UNIQUE-ID', 'UNMODIFIED-FORM'
    """
    global ECOCYC_REGULONS

    filename = os.path.join(ECOCYC_DATA_PATH, 'regulons.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_REGULONS = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'COMMON-NAME', value, True)

        set_value(entry, 'REGULATES', value)

        ECOCYC_REGULONS.append(value)

def load_protligandcplxes():
    """'CITATIONS', 'COEFFICIENT', 'COMMENT', 'COMMON-NAME', 'COMPONENT-OF',
    'COMPONENTS', 'CONSENSUS-SEQUENCE', 'CREDITS', 'DBLINKS',
    'DNA-FOOTPRINT-SIZE', 'GO-TERMS', 'LOCATIONS', 'MOLECULAR-WEIGHT-KD',
    'MOLECULAR-WEIGHT-SEQ', 'REGULATES', 'SPECIES', 'SYMMETRY', 'SYNONYMS',
    'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_PROTLIGANDCPLXES

    filename = os.path.join(ECOCYC_DATA_PATH, 'protligandcplxes.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_PROTLIGANDCPLXES = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'COMMON-NAME', value, True)

        ECOCYC_PROTLIGANDCPLXES.append(value)


def load_pathways():
    """'CITATIONS', 'CLASS-INSTANCE-LINKS', 'COMMENT', 'COMMENT-INTERNAL',
    'COMMON-NAME', 'CREDITS', 'DBLINKS', 'DIAGRAM-INFO', 'ENZYMES-NOT-USED',
    'HYPOTHETICAL-REACTIONS', 'IN-PATHWAY', 'KEY-REACTIONS',
    'PATHWAY-LINKS', 'POLYMERIZATION-LINKS', 'PREDECESSORS', 'PRIMARIES',
    'REACTION-LAYOUT', 'REACTION-LIST', 'SUB-PATHWAYS', 'SUPER-PATHWAYS',
    'SYNONYMS', 'TEMPLATE-FILE', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_PATHWAYS

    filename = os.path.join(ECOCYC_DATA_PATH, 'pathways.dat')
    data = read_ecocyc_file(filename)

    ECOCYC_PATHWAYS
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'COMMON-NAME', value, True)
        set_value(entry, 'SYNONYMS', value)

        set_value(entry, 'REACTION-LIST', value)

        set_value(entry, 'IN-PATHWAY', value)
        set_value(entry, 'SUB-PATHWAYS', value)
        set_value(entry, 'SUPER-PATHWAYS', value)

        ECOCYC_PATHWAYS.append(value)


def load_protein_features():
    """'ALTERNATE-SEQUENCE', 'ATTACHED-GROUP', 'CATALYTIC-ACTIVITY',
    'CITATIONS', 'COMMENT', 'COMMON-NAME', 'CREDITS', 'DATA-SOURCE',
    'FEATURE-OF', 'HOMOLOGY-MOTIF', 'LEFT-END-POSITION',
    'POSSIBLE-FEATURE-STATES', 'RESIDUE-NUMBER', 'RESIDUE-TYPE',
    'RIGHT-END-POSITION', 'TYPES', 'UNIQUE-ID'
    """
    global ECOCYC_PROTEIN_FEATURES

    filename = os.path.join(ECOCYC_DATA_PATH, 'protein-features.dat')
    data = read_ecocyc_file(filename)

    possible_states = []

    ECOCYC_PROTEIN_FEATURES = []
    for unique_id, entry in data.items():
        value = {}

        assert(set_value(entry, 'UNIQUE-ID', value, True))
        set_value(entry, 'TYPES', value)
        set_value(entry, 'SYNONYMS', value)
        set_value(entry, 'COMMON-NAME', value, True)

        set_value(entry, 'RESIDUE-TYPE', value, True)
        if set_value(entry, 'RESIDUE-NUMBER', value):
            value['RESIDUE-NUMBER'] = tuple(
                [int(v) for v in value['RESIDUE-NUMBER']])

        # XXX: '[LEFT,RIGHT]-END-POSITION' can take a value named 'NIL'.
        # if set_value(entry, 'LEFT-END-POSITION', value, True):
        #     value['LEFT-END-POSITION'] = int(value['LEFT-END-POSITION'])
        # if set_value(entry, 'RIGHT-END-POSITION', value, True):
        #     value['RIGHT-END-POSITION'] = int(value['RIGHT-END-POSITION'])

        set_value(entry, 'FEATURE-OF', value)
        if set_value(entry, 'POSSIBLE-FEATURE-STATES', value):
            assert(len(value['POSSIBLE-FEATURE-STATES']) == 2)
            possible_states.append(value['POSSIBLE-FEATURE-STATES'])

        ECOCYC_PROTEIN_FEATURES.append(value)

    logger.debug(set(possible_states))


# def load_cmr_gene_lookup(encoding='ISO-8859-1'):
#     def parse_braces(line, start=0):
#         tmp, retval, i = "", [], start
#         while i < len(line):
#             if line[i] == "(":
#                 if tmp != "":
#                     retval.append(tmp)
#                 i += 1
#                 tmp, i = parse_braces(line, i)
#                 retval.append(tmp)
#                 tmp = ""
#             elif line[i] == ")":
#                 break
#             elif line[i] == " ":
#                 if tmp != "":
#                     retval.append(tmp)
#                 tmp = ""
#             else:
#                 tmp += line[i]
#             i += 1
# 
#         tmp = tmp.strip()
#         if tmp != "":
#             retval.append(tmp)
#         return retval, i
# 
#     filename = os.path.join(ECOCYC_DATA_PATH, 'cmr-gene-lookup.dat')
#     with open(filename, 'r', encoding=encoding) as fin:
#         lines = fin.readlines()
#         lines = ''.join([line.strip() for line in lines])
#         data, _ = parse_braces(lines)
#         data = data[0][3]
# 
#     for key, value in data:
#         ECOCYC_CMR_GENE_LOOKUP.append((key, value.strip('"')))

# def load_gene_links(encoding='ISO-8859-1'):
#     filename = os.path.join(ECOCYC_DATA_PATH, 'gene-links.dat')
#     with open(filename, 'r', encoding=encoding) as fin:
#         reader = csv.DictReader(filter(lambda row: row[0] != '#', fin),
#                                 fieldnames=('GENE-ID', 'EG#', 'b#', 'y-name',
#                                             'CGSC-ID', 'UniProt-ID', 'GENE-NAME'),
#                                 delimiter='\t')
#         for row in reader:
#             ECOCYC_GENE_LINKS.append(row)


def load(path=None, version=None):
    global ECOCYC_DATA_PATH, ECOCYC_VERSION
    if version is not None:
        ECOCYC_VERSION = version
    if path is not None:
        ECOCYC_DATA_PATH = os.path.join(path, '{}/data'.format(ECOCYC_VERSION))

    load_classes()
    load_proteins()
    load_compounds()
    load_rnas()
    load_genes()
    load_transunits()
    load_terminators()
    load_promoters()
    load_dnabindsites()

    load_reactions()
    load_enzrxns()
    load_bindrxns()

    load_regulations()
    load_regulons()
    load_protligandcplxes()
    load_pathways()
    load_protein_features()

    # load_cmr_gene_lookup()
    # load_gene_links()

    validate_transunits()
    validate_genes()


if __name__ == '__main__':
    load()

    # print(ECOCYC_CLASSES)
    # print(ECOCYC_PROTEINS)
    # print(ECOCYC_COMPOUNDS)
    # print(ECOCYC_RNAS)
    # print(ECOCYC_GENES)
    # print(ECOCYC_TRANSUNITS)
    # print(ECOCYC_TERMINATORS)
    # print(ECOCYC_PROMOTERS)
    # print(ECOCYC_DNABINDSITES)

    all_entries = tuple(itertools.chain(
        ECOCYC_CLASSES, ECOCYC_PROTEINS, ECOCYC_COMPOUNDS, ECOCYC_RNAS,
        ECOCYC_GENES, ECOCYC_TRANSUNITS, ECOCYC_TERMINATORS, ECOCYC_PROMOTERS,
        ECOCYC_DNABINDSITES))

    genes = [entry['UNIQUE-ID'] for entry in ECOCYC_GENES]
    proteins = [entry['UNIQUE-ID'] for entry in ECOCYC_PROTEINS]
    compounds = [entry['UNIQUE-ID'] for entry in ECOCYC_COMPOUNDS]
    rnas = [entry['UNIQUE-ID'] for entry in ECOCYC_RNAS]
    class_ids = [entry['UNIQUE-ID'] for entry in ECOCYC_CLASSES]
    transunits = [entry['UNIQUE-ID'] for entry in ECOCYC_TRANSUNITS]
    terminators = [entry['UNIQUE-ID'] for entry in ECOCYC_TERMINATORS]
    unique_ids = [entry['UNIQUE-ID'] for entry in all_entries]

    bindsites = [entry['UNIQUE-ID'] for entry in itertools.chain(
        ECOCYC_DNABINDSITES, ECOCYC_PROMOTERS)]

    reactions = [entry['UNIQUE-ID'] for entry in ECOCYC_REACTIONS]
    enzrxns = [entry['UNIQUE-ID'] for entry in ECOCYC_ENZRXNS]
    regulations = [entry['UNIQUE-ID'] for entry in ECOCYC_REGULATIONS]
    pathways = [entry['UNIQUE-ID'] for entry in ECOCYC_PATHWAYS]
    features = [entry['UNIQUE-ID'] for entry in ECOCYC_PROTEIN_FEATURES]

    undefined = []

    # for entry in all_entries:
    #     if not has_key(entry, 'TYPES'):
    #         continue
    #     for t in get_value(entry, 'TYPES'):
    #         if t not in class_ids:
    #             raise ValueError('unknown type was given: ' + repr(t))

    # for entry in all_entries:
    #     if 'COMPONENTS' in entry.keys():
    #         for target in entry['COMPONENTS']:
    #             if not isinstance(target, str):
    #                 target, coef = target
    #             if (target not in unique_ids
    #                     and unwrap(target) not in unique_ids):
    #                 undefined.append(target)

    #     if 'UNMODIFIED-FORM' in entry.keys():
    #         target = entry['UNMODIFIED-FORM']
    #         if target not in unique_ids and unwrap(target) not in unique_ids:
    #             undefined.append(target)

    #     if 'GENE' in entry.keys():
    #         target = entry['GENE']
    #         if target not in unique_ids:
    #             undefined.append(target)

    #     if 'REGULATES' in entry.keys():
    #         for target in entry['REGULATES']:
    #             assert(target in regulations)

    #     if 'FEATURES' in entry.keys():
    #         for target in entry['FEATURES']:
    #             assert(target in features)

    # for entry in ECOCYC_PROMOTERS:
    #    if 'BINDS-SIGMA-FACTOR' in entry.keys():
    #        for target in entry['BINDS-SIGMA-FACTOR']:
    #            if not target in unique_ids:
    #                undefined.append(target)

    # print(sorted(set(undefined)))
    # # assert(len(undefined) == 0)

    # print(ECOCYC_REACTIONS)
    # print(ECOCYC_ENZRXNS)
    # print(ECOCYC_BINDRXNS)

    # for entry in ECOCYC_REACTIONS:
    #     for name in ('LEFT', 'RIGHT'):
    #         if name in entry.keys():
    #             for v in entry[name]:
    #                 target = v['UNIQUE-ID']
    #                 if target in itertools.chain(compounds, proteins):
    #                     continue
    #                 elif unwrap(target) in itertools.chain(
    #                         compounds, rnas, class_ids):
    #                     continue
    #                 elif target == 'E-':
    #                     continue  # XXX: exception
    #                 undefined.append(target)

    #     if 'ENZYMATIC-REACTION' in entry.keys():
    #         for target in entry['ENZYMATIC-REACTION']:
    #             if target not in enzrxns:
    #                 undefined.append(target)

    #     if 'IN-PATHWAY' in entry.keys():
    #         for target in entry['IN-PATHWAY']:
    #             if target not in itertools.chain(pathways, reactions):
    #                 undefined.append(target)
    #                 # raise ValueError(
    #                 #     "'{0}' not in pathways: ".format(target)
    #                 #     + repr(entry))

    # for entry in ECOCYC_ENZRXNS:
    #     if 'ENZYME' in entry.keys():
    #         if not entry['ENZYME'] in proteins:
    #             undefined.append(entry['ENZYME'])

    #     if 'REQUIRED-PROTEIN-COMPLEX' in entry.keys():
    #         for target in entry['REQUIRED-PROTEIN-COMPLEX']:
    #             if target not in proteins:
    #                 undefined.append(target)

    #     if 'COFACTORS' in entry.keys():
    #         for target in entry['COFACTORS']:
    #             if target in compounds:
    #                 continue
    #             elif unwrap(target) in compounds:
    #                 continue  # XXX: |Pi|
    #             elif unwrap(target) in class_ids:
    #                 continue  # XXX: |Cytochromes-B556|, NADH-P-OR-NOP
    #             elif target in proteins:
    #                 continue  # XXX: G6518-MONOMER
    #             undefined.append(target)

    #     if 'REACTION' in entry.keys():
    #         for target in entry['REACTION']:
    #             if target not in reactions:
    #                 undefined.append(target)

    # for entry in ECOCYC_BINDRXNS:
    #     if 'REACTANTS' in entry.keys():
    #         for target in entry['REACTANTS']:
    #             if target not in itertools.chain(proteins, bindsites):
    #                 # undefined.append(target)
    #                 continue  # XXX: Too many undefined ids

    #     for name in ('ACTIVATORS', 'INHIBITORS'):
    #         if name in entry.keys():
    #             for target in entry[name]:
    #                 if target[0] not in proteins:
    #                     # undefined.append(target[0])
    #                     pass # XXX: Too many undefined ids
    #                 if target[1] not in bindsites:
    #                     # undefined.append(target[1])
    #                     pass  # XXX: Too many undefined ids

    # print(sorted(set(undefined)))

    # print(ECOCYC_REGULATIONS)
    # print(ECOCYC_REGULONS)

    # for entry in ECOCYC_REGULATIONS:
    #     if 'REGULATOR' in entry.keys():
    #         for target in entry['REGULATOR']:
    #             if target in itertools.chain(proteins, rnas):
    #                 continue
    #             elif unwrap(target) in itertools.chain(compounds, class_ids):
    #                 continue
    #             pass  # XXX: Too many undefined ids
    #             # undefined.append(target)

    #     if 'REGULATED-ENTITY' in entry.keys():
    #         target = entry['REGULATED-ENTITY']
    #         if target in itertools.chain(
    #                 bindsites, enzrxns, transunits, terminators, reactions,
    #                 genes, proteins):
    #             continue
    #         undefined.append(target)

    #     # if 'ASSOCIATED-BINDING-SITE' in entry.keys():
    #     #     target = entry['ASSOCIATED-BINDING-SITE']
    #     #     if target in bindsites:
    #     #         continue
    #     #     undefined.append(target)

    #     # if 'ACCESSORY-PROTEINS' in entry.keys():
    #     #     target = entry['ACCESSORY-PROTEINS']
    #     #     if target not in proteins:
    #     #         continue
    #     #     undefined.append(target)

    # for entry in ECOCYC_PROTLIGANDCPLXES:
    #     assert(entry['UNIQUE-ID'] in proteins)

    # for entry in ECOCYC_REGULONS:
    #     assert(entry['UNIQUE-ID'] in proteins)

    #     if 'REGULATES' in entry.keys():
    #         for target in entry['REGULATES']:
    #             assert(target in regulations)

    # print(ECOCYC_PATHWAYS)

    # for entry in ECOCYC_PATHWAYS:
    #     for name in ('SUB-PATHWAYS', 'SUPER-PATHWAYS', 'IN-PATHWAY'):
    #         if name in entry.keys():
    #             for target in entry[name]:
    #                 assert(target in pathways)
    #     if 'REACTION-LIST' in entry.keys():
    #         for target in entry['REACTION-LIST']:
    #             if target not in itertools.chain(reactions, pathways):
    #                 undefined.append(target)

    print(ECOCYC_PROTEIN_FEATURES[: 10])

    for entry in ECOCYC_PROTEIN_FEATURES:
        if 'FEATURE-OF' in entry.keys():
            for target in entry['FEATURE-OF']:
                if target not in proteins:
                    undefined.append(target)

    undefined = sorted(set(undefined))
    print(undefined, len(undefined))
