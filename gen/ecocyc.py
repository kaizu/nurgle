#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import csv
# import copy
# import os.path
# import itertools
import logging
log_ = logging.getLogger(__name__)
from enum import Enum, auto
# from . import utils

from . import _ecocyc

__program__ = 'ecocyc.py'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu'


ECOCYC_LOADED = False

def load(*args, **kwargs):
    global ECOCYC_LOADED
    if not ECOCYC_LOADED:
        _ecocyc.load(*args, **kwargs)
        ECOCYC_LOADED = True

def find_entry(root, entry_id):
    entries = []
    for r in root:
        if r['UNIQUE-ID'] == entry_id:
            entries.append(r)
    if len(entries) == 0:
        # log_.warn('[{}] could not be found.'.format(entry_id))
        return None
    elif len(entries) > 1:
        # log_.warn('[{}] has multiple entries.'.format(entry_id))
        return None
    return entries[0]

def find_reaction(entry_id):
    return find_entry(_ecocyc.ECOCYC_REACTIONS, entry_id)

def find_enzrxn(entry_id):
    return find_entry(_ecocyc.ECOCYC_ENZRXNS, entry_id)

def proteins():
    return _ecocyc.ECOCYC_PROTEINS

def rnas():
    return _ecocyc.ECOCYC_RNAS

def transunits():
    return _ecocyc.ECOCYC_TRANSUNITS

class EcocycKind(Enum):
    UNKNOWN = auto()
    REACTION = auto()
    ENZRXN = auto()
    PROTEIN = auto()
    RNA = auto()
    TRANSUNIT = auto()

def kind(entry_id):
    for root, res in (
        (_ecocyc.ECOCYC_REACTIONS, EcocycKind.REACTION),
        (_ecocyc.ECOCYC_ENZRXNS, EcocycKind.ENZRXN),
        (_ecocyc.ECOCYC_PROTEINS, EcocycKind.PROTEIN),
        (_ecocyc.ECOCYC_RNAS, EcocycKind.RNA),
        (_ecocyc.ECOCYC_TRANSUNITS, EcocycKind.TRANSUNIT),
        ):
        if find_entry(root, entry_id) is not None:
            return res
    return EcocycKind.UNKNOWN