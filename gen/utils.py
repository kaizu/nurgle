#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path
import logging

try:
    from urllib.request import Request, urlopen, URLError  # Python3
except ImportError:
    from urllib2 import Request, urlopen, URLError  # Python2


__program__ = "utils.py"
__version__ = "0.1"
__author__ = "Kazunari Kaizu"

# DATA_PATH = os.path.abspath(os.path.join(
#     os.path.dirname(os.path.abspath(__file__)), '../data'))

__all__ = ['collect_keys', 'collect_values', 'collect_if']

logger = logging.getLogger('libecoli')


def isiterable(value):
    return isinstance(value, (list, tuple, set, frozenset))


def collect_keys(data):
    retval = set()
    for elem in data:
        retval |= set(elem.keys())
    return tuple(retval)


def collect_values(data, key):
    return tuple(set([elem[key] for elem in data if key in elem.keys()]))


def collect_values_with_iterable(collection, key):
    retval = []
    for doc in collection:
        value = doc.get(key)
        if value is None:
            continue
        elif isiterable(value):
            retval.extend(value)
        else:
            retval.append(value)
    return tuple(set(retval))


def collect_if(data, cond):
    if isinstance(cond, dict):
        def predicator(elem):
            for k, v in cond.items():
                if elem.get(k) != v:
                    return False
            return True
        return tuple([elem for elem in data if predicator(elem)])
    return tuple([elem for elem in data if cond(elem)])


def fetch_url(url, filename, encoding='utf-8', errors='strict'):
    req = Request(url)

    try:
        response = urlopen(req)
    except URLError as e:
        if hasattr(e, 'reason'):
            logger.error('We failed to reach a server.')
            logger.error('Reason: ' + e.reason)
            return False
        elif hasattr(e, 'code'):
            logger.error('The server couldn\'t fulfill the request.')
            logger.error('Error code: ' + e.code)
            return False
    else:
        # everything is fine
        filename = os.path.abspath(filename)
        root_path = os.path.dirname(filename)
        if not os.path.isdir(root_path):
            os.makedirs(root_path)
            logger.info('Made directories [{0}]'.format(root_path))

        logger.info('Dowloading [{0}] from [{1}] ...'.format(filename, url))
        if encoding != 'binary':
            with open(filename, 'w') as fout:
                fout.write(response.read().decode(encoding, errors))
        else:
            with open(filename, 'wb') as fout:
                fout.write(response.read())
    return True


def indexing(collection):
    retval = collections.defaultdict(list)
    for i, doc in enumerate(collection):
        for value in doc.values():
            if isiterable(value):
                for value_ in value:
                    retval[value_].append(i)
            else:
                retval[value].append(i)
    for key in retval.keys():
        retval[key] = tuple(set(retval[key]))
    return retval


def find_any(collection, key, value=None, collection_index=None):
    if collection_index is None:
        retval = []
        if value is None:
            for doc in collection:
                for value_ in doc.values():
                    if key == value_ or (
                        isiterable(value_) and key in value_):
                        retval.append(doc)
        else:
            for doc in collection:
                value_ = doc.get(key)
                if value == value_ or (
                    isiterable(value_) and value in value_):
                    retval.append(doc)
        return retval
    else:
        if value is None:
            indices = collection_index.get(key)
            return ([collection[i] for i in indices]
                if indices is not None else [])
        else:
            indices = collection_index.get(value)
            if indices is None:
                return []
            retval = []
            for i in indices:
                value_ = collection[i].get(key)
                if value == value_ or (
                    isiterable(value_) and value in value_):
                    retval.append(doc)
            return retval
