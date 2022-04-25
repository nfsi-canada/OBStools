from obstools.atacr import DayNoise, StaNoise, TFNoise, EventStream
from obstools.comply import Comply
from . import get_meta
import pytest


def test_daynoise_demo():
    return DayNoise('demo')


def test_day_QC(tmp_path):
    daynoise = test_daynoise_demo()
    daynoise.QC_daily_spectra(fig_QC=True, save=tmp_path)
    return daynoise


def test_day_ncomp_opts(tmp_path):
    daynoise = test_daynoise_demo()
    daynoise.ncomp = 3
    daynoise.QC_daily_spectra(fig_QC=True, save=tmp_path)
    daynoise.ncomp = 2
    daynoise.QC_daily_spectra(fig_QC=True, save=tmp_path, smooth=False)


def test_day_average(tmp_path):
    daynoise = test_daynoise_demo()
    daynoise.average_daily_spectra(
        fig_average=True, fig_coh_ph=True,
        save=tmp_path)
    return daynoise


def test_day_save(tmp_path):
    daynoise = test_daynoise_demo()
    d = tmp_path / "tmp"
    daynoise.save(d)


def test_stanoise_demo():
    return StaNoise('demo')


def test_stanoise_day_demo():
    daynoise = test_daynoise_demo()
    stanoise = StaNoise(daylist=[daynoise])
    return stanoise


def test_sta_QC(tmp_path):
    stanoise = test_stanoise_demo()
    stanoise.QC_sta_spectra(fig_QC=True, save=tmp_path)
    return stanoise


def test_sta_ncomp_opts(tmp_path):
    stanoise = test_stanoise_demo()
    stanoise.ncomp = 3
    stanoise.QC_sta_spectra(fig_QC=True, save=tmp_path)
    stanoise = test_stanoise_demo()
    stanoise.ncomp = 2
    stanoise.initialized = None
    stanoise.QC_sta_spectra(fig_QC=True, save=tmp_path)


def test_sta_average(tmp_path):
    stanoise = test_stanoise_demo()
    with pytest.raises(Exception):
        assert stanoise.QC_sta_spectra()

    stanoise.average_sta_spectra(fig_average=True, save=tmp_path)
    return stanoise


def test_sta_operations():
    dn1 = test_daynoise_demo()
    dn2 = test_daynoise_demo()
    sn = StaNoise()
    sn += dn1
    sn.append(dn2)
    sn.extend([dn1, dn2])
    with pytest.raises(TypeError):
        assert sn.append([dn1, dn2])
    sn.extend(sn)
    with pytest.raises(TypeError):
        assert sn.extend([float])
    with pytest.raises(TypeError):
        assert sn.extend({})
    sn = StaNoise(daylist=[dn1])
    with pytest.raises(Exception):
        assert sn.init()


def test_sta_save(tmp_path):
    stanoise = test_stanoise_demo()
    d = tmp_path / "tmp"
    with pytest.raises(Exception):
        assert stanoise.save(d)

    stanoise.average_sta_spectra()
    stanoise.save(d)


def test_evstream_demo(tmp_path):
    eventstream = EventStream('demo')
    d = tmp_path / "tmp"
    eventstream.save(d)
    return eventstream


def test_tfnoise_day_demo(tmp_path):
    daynoise = test_daynoise_demo()
    with pytest.raises(TypeError):
        assert TFNoise([])
    with pytest.raises(Exception):
        TFNoise(objnoise=daynoise)
    daynoise.average_daily_spectra()
    tfnoise_day = TFNoise(daynoise)
    tfnoise_day.transfer_func()
    d = tmp_path / "tmp"
    tfnoise_day.save(d)
    return tfnoise_day


def test_tfnoise_sta_demo(tmp_path):
    stanoise = test_sta_average(tmp_path)
    tfnoise_sta = TFNoise(stanoise)
    tfnoise_sta.transfer_func()
    return tfnoise_sta


def test_comply_day_demo(tmp_path):
    daynoise = test_day_average(tmp_path)
    sta = get_meta.get_stdb()
    comply_day = Comply(objnoise=daynoise, elev=sta.elevation*1.e3)
    comply_day.calculate_compliance()
    return comply_day


def test_comply_sta_demo(tmp_path):
    stanoise = test_sta_average(tmp_path)
    sta = get_meta.get_stdb()
    comply_sta = Comply(objnoise=stanoise, elev=sta.elevation*1.e3)
    comply_sta.calculate_compliance()
    return comply_sta


def test_comply_fail(tmp_path):
    import pytest
    import pickle
    import os

    with pytest.raises(Exception):
        assert Comply()

    sta = get_meta.get_stdb()
    with pytest.raises(Exception):
        assert Comply(elev=sta.elevation*1.e3, objnoise=[])

    objnoise = test_day_average(tmp_path)
    objnoise.av = None
    with pytest.raises(Exception):
        assert Comply(elev=sta.elevation*1.e3, objnoise=objnoise)

    daynoise = test_day_average(tmp_path)
    comply_day = Comply(objnoise=daynoise, elev=sta.elevation*1.e3)
    comply_day.calculate_compliance()
    d = tmp_path / "tmp"
    comply_day.save(d, form='pkl')
    dd = d.parent / (d.name + '.' + 'pkl')
    assert dd.exists()

    d = tmp_path / "tmp"
    comply_day.save(d, form='csv')
    dd = d.parent / (d.name + '.' + 'csv')
    assert dd.exists()

    del comply_day.complyfunc['ZP-H']

    d = tmp_path / "tmp"
    comply_day.save(d, form='csv')
    dd = d.parent / (d.name + '.' + 'csv')
    assert dd.exists()

    del comply_day.complyfunc['ZP-21']

    d = tmp_path / "tmp"
    comply_day.save(d, form='csv')
    dd = d.parent / (d.name + '.' + 'csv')
    assert dd.exists()
