def test_stdb_import():
    import stdb

def test_obspy_import():
    import obspy

def test_obstools_modules():
    import obstools
    from obstools.atacr import classes
    from obstools.atacr import utils
    import matplotlib
    matplotlib.use('Agg')
    from obstools.atacr import plot
    from obstools.atacr.classes import DayNoise, StaNoise, TFNoise, EventStream, Power, Cross, Rotation 