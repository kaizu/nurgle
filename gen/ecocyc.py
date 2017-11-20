#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import csv
# import copy
# import os.path
# import itertools
import logging
log_ = logging.getLogger(__name__)

# from . import utils

from . import _ecocyc

__program__ = 'ecocyc.py'
__version__ = '0.1'
__author__ = 'Kazunari Kaizu'


load = _ecocyc.load

def find_entry(root, entry_id):
    entries = []
    for r in root:
        if r['UNIQUE-ID'] == entry_id:
            entries.append(r)
    if len(entries) == 0:
        log_.warn('[{}] could not be found.'.format(entry_id))
        return None
    elif len(entries) > 1:
        log_.warn('[{}] has multiple entries.'.format(entry_id))
        return None
    return entries[0]

def find_reaction(entry_id):
    return find_entry(_ecocyc.ECOCYC_REACTIONS, entry_id)

def find_enzrxn(entry_id):
    return find_entry(_ecocyc.ECOCYC_ENZRXNS, entry_id)
