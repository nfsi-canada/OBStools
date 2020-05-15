from obspy import UTCDateTime
from obstools.atacr import utils
from pkg_resources import resource_filename
from pathlib import Path

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

def test_get_event():
    eventpath = exmpl_path / 'event'
    tstart = UTCDateTime('2012-04-09')
    tend = UTCDateTime('2012-04-10')
    tr1, tr2, trZ, trP = utils.get_event(eventpath, tstart, tend)

