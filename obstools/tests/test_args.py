import numpy as np
from obstools.atacr import arguments
from pkg_resources import resource_filename
from pathlib import Path


dbfile = resource_filename('obstools',
                           'examples/meta/M08A.pkl')

def test_get_daylong_arguments():
    args = arguments.get_daylong_arguments([dbfile])
    return args

def test_get_event_arguments():
    args = arguments.get_event_arguments([dbfile])
    return args

def test_get_dailyspec_arguments():
    args = arguments.get_dailyspec_arguments([dbfile])
    return args

def test_get_cleanspec_arguments():
    args = arguments.get_cleanspec_arguments([dbfile])
    return args

def test_get_transfer_arguments():
    args = arguments.get_transfer_arguments([dbfile])
    return args

def test_get_correct_arguments():
    args = arguments.get_correct_arguments([dbfile])
    return args

# def test_dirs(tmp_path):
#     args = test_dl_calc_args()
#     db = get_meta.get_stdb()
#     stkey = 'YH.LOBS3'

#     outdir = tmp_path / args.saveloc 
#     if not outdir.exists():
#         outdir.mkdir()

#     outdir = outdir / stkey.upper()

#     if not outdir.exists():
#         outdir.mkdir()

