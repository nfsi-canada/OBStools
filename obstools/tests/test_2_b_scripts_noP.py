import stdb
import numpy as np
import shutil
from pathlib import Path
from pkg_resources import resource_filename
from obspy.clients.fdsn import Client
from obstools.atacr import DayNoise, StaNoise, TFNoise
from obstools.atacr import EventStream, Power, Cross, Rotation
from obstools.atacr import utils, plotting

dbfile = resource_filename('obstools',
                           'examples/meta/M08A.pkl')

curdir = Path.cwd()
datadir = curdir / 'DATA'
avgdir = curdir / 'AVG_STA'
cmpdir = curdir / 'COMPL_STA'
evdir = curdir / 'EVENTS'
specdir = curdir / 'SPECTRA'
tfdir = curdir / 'TF_STA'


# Test with horizontals only (no P)
def test_11_data_noP():
    from obstools.scripts import atacr_download_data as atacr
    args0 = atacr.get_daylong_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--sampling-rate', '1.0', '--channels', '12'])
    atacr.main(args=args0)


def test_12_daily_noP():
    from obstools.scripts import atacr_daily_spectra as atacr
    args0 = atacr.get_dailyspec_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--figQC',
        '--figAverage', '--save-fig'])
    atacr.main(args=args0)


def test_13_clean_noP():
    from obstools.scripts import atacr_clean_spectra as atacr
    args0 = atacr.get_cleanspec_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--save-fig',
        '--figCross', '--figCoh'])
    atacr.main(args=args0)


def test_14_trans_noP():
    from obstools.scripts import atacr_transfer_functions as atacr
    args0 = atacr.get_transfer_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--save-fig', '--figTF'])
    atacr.main(args=args0)


def test_15_event_noP():
    from obstools.scripts import atacr_download_event as atacr
    args0 = atacr.get_event_arguments([
        dbfile, '--keys', '7D.M08A', '-O',
        '--start', '2012-03-08', '--end', '2012-03-10',
        '--min-mag', '6.3', '--max-mag', '6.7', '--window', '7200.',
        '--sampling-rate', '1.0', '--channels', '12'])
    atacr.main(args=args0)


def test_16_correct_noP():
    from obstools.scripts import atacr_correct_event as atacr
    args0 = atacr.get_correct_arguments([
        dbfile, '--keys', '7D.M08A', '-O', '--figRaw',
        '--figClean', '--save-fig', '--save'])
    atacr.main(args=args0)


def test_17_rmtree():
    shutil.rmtree(datadir)
    shutil.rmtree(avgdir)
    shutil.rmtree(evdir)
    shutil.rmtree(specdir)
    shutil.rmtree(tfdir)
