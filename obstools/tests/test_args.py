import numpy as np
from pkg_resources import resource_filename
from pathlib import Path
import pytest


dbfile = resource_filename('obstools',
                           'examples/meta/M08A.pkl')


def test_get_daylong_arguments():
    from obstools.scripts import atacr_download_data as atacr
    # no stdb
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments()
    # defaults
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.MM08',
        '--start', '2012-03-01', '--end', '2012-03-05'])
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--channels', 'J,P'])
    # start time
    args = atacr.get_daylong_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_daylong_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--end', 'abcd'])
    # user auth
    args = atacr.get_daylong_arguments([
        dbfile, '-U', 'user:name'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '-U', 'abcd'])
    # sampling rate
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--sampling-rate', 'abcd'])
    # units
    with pytest.raises(Exception):
        assert atacr.get_daylong_arguments([
            dbfile, '--units', 'abcd'])
    # pre-filt
    with pytest.raises(Exception):
        assert atacr.get_daylong_arguments([
            dbfile, '--pre-filt', '0.1,0.2,0.3'])

    return args0


def test_get_event_arguments():
    from obstools.scripts import atacr_download_event as atacr
    # defaults
    args = atacr.get_event_arguments([dbfile])
    # no stdb
    with pytest.raises(SystemExit):
        assert atacr.get_event_arguments()
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A', '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--keys', 'MM08,M08A', '--channels', 'J,P'])
    return args


def test_get_dailyspec_arguments():
    from obstools.scripts import atacr_daily_spectra as atacr
    args = atacr.get_dailyspec_arguments([dbfile])
    # no stdb
    with pytest.raises(SystemExit):
        assert atacr.get_dailyspec_arguments()
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A', '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--keys', 'MM08,M08A', '--channels', 'J,P'])
    return args


def test_get_cleanspec_arguments():
    from obstools.scripts import atacr_clean_spectra as atacr
    args = atacr.get_cleanspec_arguments([dbfile])
    # no stdb
    with pytest.raises(SystemExit):
        assert atacr.get_cleanspec_arguments()
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A', '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--keys', 'MM08,M08A', '--channels', 'J,P'])
    return args


def test_get_transfer_arguments():
    from obstools.scripts import atacr_transfer_functions as atacr
    args = atacr.get_transfer_arguments([dbfile])
    # no stdb
    with pytest.raises(SystemExit):
        assert atacr.get_transfer_arguments()
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A', '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--keys', 'MM08,M08A', '--channels', 'J,P'])
    return args


def test_get_correct_arguments():
    from obstools.scripts import atacr_correct_event as atacr
    args = atacr.get_correct_arguments([dbfile])
    # no stdb
    with pytest.raises(SystemExit):
        assert atacr.get_correct_arguments()
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A', '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--keys', 'MM08,M08A', '--channels', 'J,P'])
    return args


def test_get_comply():
    from obstools.scripts import comply_calculate as comply
    args = comply.get_comply_arguments([dbfile])
    # no stdb
    with pytest.raises(SystemExit):
        assert comply.get_comply_arguments()
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08,M08A', '--channels', 'H,P'])
    with pytest.raises(SystemExit):
        assert atacr.get_daylong_arguments([
            dbfile, '--keys', 'MM08,M08A', '--channels', 'J,P'])
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
