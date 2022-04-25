import numpy as np
from pkg_resources import resource_filename
from pathlib import Path
import pytest
from parser import ParserError


dbfile = resource_filename('obstools',
                           'examples/meta/M08A.pkl')


def test_get_daylong_arguments():
    from obstools.scripts import atacr_download_data as atacr
    # no stdb
    with pytest.raises(SystemExit):
        atacr.get_daylong_arguments([])
    # defaults
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.MM08',
        '--start', '2012-03-01', '--end', '2012-03-05'])
    # keys
    args = atacr.get_daylong_arguments([
        dbfile, '--keys', 'MM08'])
    # channels
    args = atacr.get_daylong_arguments([
        dbfile, '--channels', '12,P'])
    with pytest.raises(SystemExit):
        atacr.get_daylong_arguments([
            dbfile, '--channels', 'J,P'])
    # start time
    args = atacr.get_daylong_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_daylong_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_daylong_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_daylong_arguments([
            dbfile, '--end', 'abcd'])
    # user auth
    args = atacr.get_daylong_arguments([
        dbfile, '-U', 'user:name'])
    with pytest.raises(SystemExit):
        atacr.get_daylong_arguments([
            dbfile, '-U', 'abcd'])
    # sampling rate
    with pytest.raises(SystemExit):
        atacr.get_daylong_arguments([
            dbfile, '--sampling-rate', 'abcd'])
    # units
    with pytest.raises(Exception):
        atacr.get_daylong_arguments([
            dbfile, '--units', 'abcd'])
    # pre-filt
    with pytest.raises(Exception):
        atacr.get_daylong_arguments([
            dbfile, '--pre-filt', '0.1,0.2,0.3'])

    return args0


def test_get_event_arguments():
    from obstools.scripts import atacr_download_event as atacr
    # no stdb
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([])
    # defaults
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.MM08',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7'])
    # keys
    args = atacr.get_event_arguments([
        dbfile, '--keys', '7D.MM08'])
    # channels
    args = atacr.get_event_arguments([
        dbfile, '--channels', '12,P'])
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--channels', 'J,P'])
    # start time
    args = atacr.get_event_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_event_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--end', 'abcd'])
    # user auth
    args = atacr.get_event_arguments([
        dbfile, '-U', 'user:name'])
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '-U', 'abcd'])
    # sampling rate
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--sampling-rate', 'abcd'])
    # magnitudes
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--min-mag', 'abcd'])
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--max-mag', 'abcd'])
    # distances
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--min-dist', 'abcd'])
    with pytest.raises(SystemExit):
        atacr.get_event_arguments([
            dbfile, '--max-dist', 'abcd'])
    # pre-filt
    with pytest.raises(Exception):
        atacr.get_event_arguments([
            dbfile, '--pre-filt', '0.1,0.2,0.3'])
    return args0


def test_get_dailyspec_arguments():
    from obstools.scripts import atacr_daily_spectra as atacr
    # no stdb
    with pytest.raises(SystemExit):
        atacr.get_dailyspec_arguments([])
    # defaults
    args0 = atacr.get_dailyspec_arguments([dbfile])
    # keys
    args = atacr.get_dailyspec_arguments([
        dbfile, '--keys', '7D.MM08'])
    # raw or smooth
    args = atacr.get_dailyspec_arguments([
        dbfile, '--raw'])
    # start time
    args = atacr.get_dailyspec_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_dailyspec_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_dailyspec_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_dailyspec_arguments([
            dbfile, '--end', 'abcd'])
    # frequencies
    args = atacr.get_dailyspec_arguments([
        dbfile, '--freq-band', '0.004,2.'])
    with pytest.raises(Exception):
        atacr.get_dailyspec_arguments([
            dbfile, '--freq-band', '0.1'])

    return args0


def test_get_cleanspec_arguments():
    from obstools.scripts import atacr_clean_spectra as atacr
    # no stdb
    with pytest.raises(SystemExit):
        atacr.get_cleanspec_arguments([])
    # defaults
    args0 = atacr.get_cleanspec_arguments([dbfile])
    # keys
    args = atacr.get_cleanspec_arguments([
        dbfile, '--keys', '7D.MM08'])
    # start time
    args = atacr.get_cleanspec_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_cleanspec_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_cleanspec_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_cleanspec_arguments([
            dbfile, '--end', 'abcd'])
    # frequencies
    args = atacr.get_cleanspec_arguments([
        dbfile, '--freq-band', '0.004,2.'])
    with pytest.raises(Exception):
        atacr.get_cleanspec_arguments([
            dbfile, '--freq-band', '0.1'])
    return args0


def test_get_transfer_arguments():
    from obstools.scripts import atacr_transfer_functions as atacr
    # no stdb
    with pytest.raises(SystemExit):
        atacr.get_transfer_arguments([])
    # defaults
    args0 = atacr.get_transfer_arguments([dbfile])
    # keys
    args = atacr.get_transfer_arguments([
        dbfile, '--keys', '7D.MM08'])
    # start time
    args = atacr.get_transfer_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_transfer_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_transfer_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_transfer_arguments([
            dbfile, '--end', 'abcd'])
    # skip clean
    with pytest.raises(SystemExit):
        atacr.get_transfer_arguments([
            dbfile, '--skip-clean', '--skip-daily'])
    return args0


def test_get_correct_arguments():
    from obstools.scripts import atacr_correct_event as atacr
    # no stdb
    with pytest.raises(SystemExit):
        atacr.get_correct_arguments([])
    # defaults
    args0 = atacr.get_correct_arguments([dbfile])
    # keys
    args = atacr.get_correct_arguments([
        dbfile, '--keys', '7D.MM08'])
    # start time
    args = atacr.get_correct_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_correct_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = atacr.get_correct_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        atacr.get_correct_arguments([
            dbfile, '--end', 'abcd'])
    # skip clean
    with pytest.raises(SystemExit):
        atacr.get_correct_arguments([
            dbfile, '--skip-clean', '--skip-daily'])
    return args0


def test_get_comply_arguments():
    from obstools.scripts import comply_calculate as comply
    # no stdb
    with pytest.raises(SystemExit):
        comply.get_comply_arguments([])
    # defaults
    args0 = comply.get_comply_arguments([dbfile])
    # keys
    args = comply.get_comply_arguments([
        dbfile, '--keys', '7D.MM08'])
    # start time
    args = comply.get_comply_arguments([
        dbfile, '--start', '2020-01-01'])
    with pytest.raises(SystemExit):
        comply.get_comply_arguments([
            dbfile, '--start', 'abcd'])
    # end time
    args = comply.get_comply_arguments([
        dbfile, '--end', '2020-01-01'])
    with pytest.raises(SystemExit):
        comply.get_comply_arguments([
            dbfile, '--end', 'abcd'])
    # skip clean
    with pytest.raises(SystemExit):
        comply.get_comply_arguments([
            dbfile, '--skip-clean', '--skip-daily'])
    # save format
    with pytest.raises(SystemExit):
        comply.get_comply_arguments([
            dbfile, '--save-format', 'abcd'])
    return args0

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
