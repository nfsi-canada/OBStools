from obspy import UTCDateTime, read
from obstools.atacr import utils
from pkg_resources import resource_filename
from pathlib import Path
import shutil
import glob
import os

# Path where data are located
exmpl_path = Path(resource_filename('obstools','examples'))

def test_get_data():
    datapath = exmpl_path / 'data'
    tstart = UTCDateTime('2012-04-01')
    tend = UTCDateTime('2012-04-04')
    trN1, trN2, trNZ, trNP = utils.get_data(datapath, tstart, tend)
    assert trN1 is not None
    assert trN2 is not None
    assert trNZ is not None
    assert trNP is not None
    return trN1, trN2, trNZ, trNP

def test_get_H_data(tmp_path):
    datapath = exmpl_path / 'data'

    # Test only H data
    for filename in glob.glob(os.path.join(datapath, '*1.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*2.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-04-01')
    tend = UTCDateTime('2012-04-04')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)
    assert trN1 is not None
    assert trN2 is not None
    assert trNZ is not None
    assert len(trNP) == 0

def test_get_P_data(tmp_path):
    datapath = exmpl_path / 'data'

    # Test only P data
    for filename in glob.glob(os.path.join(datapath, '*H.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-04-01')
    tend = UTCDateTime('2012-04-04')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)
    assert len(trN1) == 0
    assert len(trN2) == 0
    assert trNZ is not None
    assert trNP is not None

    for filename in glob.glob(os.path.join(datapath, '*H.SAC')):
        stP = read(filename)
        stP[0].resample(1.)
        stP[0].write(filename, format='SAC')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)

    for filename in glob.glob(os.path.join(datapath, '*H.SAC')):
        stP = read(filename)
        stP[0].resample(10.)
        stP[0].write(filename, format='SAC')
    trN1, trN2, trNZ, trNP = utils.get_data(tmp_path, tstart, tend)

def test_get_event():
    eventpath = exmpl_path / 'event'
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(eventpath, tstart, tend)
    assert tr1 is not None
    assert tr2 is not None
    assert trZ is not None
    assert trP is not None
    return tr1, tr2, trZ, trP

def test_get_H_event(tmp_path):
    datapath = exmpl_path / 'event'

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
    assert tr1 is not None
    assert tr2 is not None
    assert trZ is not None
    assert len(trP) == 0

def test_get_P_event(tmp_path):
    datapath = exmpl_path / 'event'

    # Test only P data
    for filename in glob.glob(os.path.join(datapath, '*H.SAC')):
        shutil.copy(filename, tmp_path)
    for filename in glob.glob(os.path.join(datapath, '*Z.SAC')):
        shutil.copy(filename, tmp_path)
    tstart = UTCDateTime('2012-03-08')
    tend = UTCDateTime('2012-03-10')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)

    filename = glob.glob(os.path.join(datapath, '*H.SAC'))[0]
    stP = read(filename)
    print(stP)
    stP[0].resample(1.)
    stP[0].write(filename, format='SAC')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)

    stP[0].resample(10.)
    stP[0].write(filename, format='SAC')
    tr1, tr2, trZ, trP = utils.get_event(tmp_path, tstart, tend)

