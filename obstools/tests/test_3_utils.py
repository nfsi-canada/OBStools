from obspy import UTCDateTime, read
from obstools.atacr import utils
from pkg_resources import resource_filename
from pathlib import Path
import shutil
import glob
import os

# Path where data are located
exmpl_path = Path(resource_filename('obstools', 'examples'))


def test_get_data():
    datapath = Path('DATA') / '7D.M08A'
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-10')
    trN1, trN2, trNZ, trNP = utils.get_data(datapath, tstart, tend)
    assert trN1 is not None
    assert trN2 is not None
    assert trNZ is not None
    assert trNP is not None
    return trN1, trN2, trNZ, trNP


def test_get_H_data(tmp_path):
    datapath = Path('DATA') / '7D.M08A'

    # Test only H data
    for filename in glob.glob(os.path.join(datapath, '*1.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*2.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-10')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)
    assert len(trN1) > 0
    assert len(trN2) > 0
    assert len(trNZ) > 0
    assert [len(tr.data) == 0 for tr in trNP]


def test_get_P_data(tmp_path):
    datapath = Path('DATA') / '7D.M08A'

    # Test only P data
    for filename in glob.glob(os.path.join(datapath, '*H.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-10')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)
    assert [len(tr.data) == 0 for tr in trN1]
    assert [len(tr.data) == 0 for tr in trN2]
    assert len(trNZ) > 0
    assert len(trNP) > 0


def test_get_P_data_sr1(tmp_path):
    datapath = Path('DATA') / '7D.M08A'

    for filename in glob.glob(os.path.join(datapath, '*.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(tmp_path, '*H.SAC')):
        stP = read(filename)
        stP[0].resample(0.5)
        stP[0].write(filename, format='SAC')
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-10')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)


def test_get_P_data_sr2(tmp_path):
    datapath = Path('DATA') / '7D.M08A'

    for filename in glob.glob(os.path.join(datapath, '*.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(tmp_path, '*H.SAC')):
        stP = read(filename)
        stP[0].resample(10.)
        stP[0].write(filename, format='SAC')
    tstart = UTCDateTime('2012-03-01')
    tend = UTCDateTime('2012-03-10')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)


def test_get_event():
    datapath = Path('EVENTS') / '7D.M08A'
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(datapath, tstart, tend)
    assert len(tr1) > 0
    assert len(tr2) > 0
    assert len(trZ) > 0
    assert len(trP) > 0
    return tr1, tr2, trZ, trP


def test_get_H_event(tmp_path):
    datapath = Path('EVENTS') / '7D.M08A'

    # Test only H data
    for filename in glob.glob(os.path.join(datapath, '*1.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*2.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)
    assert len(tr1) > 0
    assert len(tr2) > 0
    assert len(trZ) > 0
    assert [len(tr.data) == 0 for tr in trP]


def test_get_P_event(tmp_path):
    datapath = Path('EVENTS') / '7D.M08A'

    # Test only P data
    for filename in glob.glob(os.path.join(datapath, '*H.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)
    assert len(trZ) > 0
    assert len(trP) > 0
    assert [len(tr.data) == 0 for tr in tr1]
    assert [len(tr.data) == 0 for tr in tr2]


def test_get_P_event_sr1(tmp_path):
    datapath = Path('EVENTS') / '7D.M08A'

    for filename in glob.glob(os.path.join(datapath, '*.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(tmp_path, '*H.SAC')):
        stP = read(filename)
        stP[0].resample(0.5)
        stP[0].write(filename, format='SAC')
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)


def test_get_P_event_sr2(tmp_path):
    datapath = Path('EVENTS') / '7D.M08A'

    for filename in glob.glob(os.path.join(datapath, '*.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(tmp_path, '*H.SAC')):
        stP = read(filename)
        stP[0].resample(10.)
        stP[0].write(filename, format='SAC')
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)
