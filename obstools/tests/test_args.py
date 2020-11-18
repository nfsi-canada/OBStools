import numpy as np
from pkg_resources import resource_filename
from pathlib import Path


dbfile = resource_filename('obstools',
                           'examples/meta/M08A.pkl')

def test_get_daylong_arguments():
    from obstools.scripts import atacr_download_data as atacr
    args = atacr.get_daylong_arguments([dbfile])
    return args

def test_get_event_arguments():
    from obstools.scripts import atacr_download_event as atacr
    args = atacr.get_event_arguments([dbfile])
    return args

def test_get_dailyspec_arguments():
    from obstools.scripts import atacr_daily_spectra as atacr
    args = atacr.get_dailyspec_arguments([dbfile])
    return args

def test_get_cleanspec_arguments():
    from obstools.scripts import atacr_clean_spectra as atacr
    args = atacr.get_cleanspec_arguments([dbfile])
    return args

def test_get_transfer_arguments():
    from obstools.scripts import atacr_transfer_functions as atacr
    args = atacr.get_transfer_arguments([dbfile])
    return args

def test_get_correct_arguments():
    from obstools.scripts import atacr_correct_event as atacr
    args = atacr.get_correct_arguments([dbfile])
    return args

def test_get_comply():
    from obstools.scripts import comply_calculate as comply
    args = comply.get_comply_arguments([dbfile])
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

