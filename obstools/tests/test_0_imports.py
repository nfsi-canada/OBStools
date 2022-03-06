def test_stdb_import():
    import stdb


def test_obspy_import():
    import obspy


def test_obstools_modules():
    import obstools
    import obstools.atacr
    import obstools.comply
    from obstools.atacr import classes
    from obstools.atacr import utils
    import matplotlib
    matplotlib.use('Agg')
    from obstools.atacr import plotting
    from obstools.atacr.classes import DayNoise, StaNoise, TFNoise
    from obstools.atacr.classes import EventStream, Power, Cross, Rotation
    from obstools.comply.classes import Comply
