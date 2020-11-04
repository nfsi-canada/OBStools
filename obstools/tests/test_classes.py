from obstools.atacr import DayNoise, StaNoise, TFNoise, EventStream
from obstools.comply import Comply
from . import get_meta

def test_daynoise_demo():
    return DayNoise('demo')

def test_day_QC():
    daynoise = test_daynoise_demo()
    daynoise.QC_daily_spectra()
    return daynoise

def test_day_average():
    daynoise = test_day_QC()
    daynoise.average_daily_spectra()
    return daynoise

def test_stanoise_demo():
    return StaNoise('demo')

def test_stanoise_day_demo():
    daynoise = test_daynoise_demo()
    stanoise = StaNoise(daylist=[daynoise])
    return stanoise

def test_sta_QC():
    stanoise = test_stanoise_demo()
    stanoise.QC_sta_spectra()
    return stanoise

def test_sta_average():
    stanoise = test_sta_QC()
    stanoise.average_sta_spectra()
    return stanoise

def test_evstream_demo():
    return EventStream('demo')

def test_tfnoise_day_demo():
    daynoise = test_day_average()
    tfnoise_day = TFNoise(daynoise)
    tfnoise_day.transfer_func()
    return tfnoise_day

def test_tfnoise_sta_demo():
    stanoise = test_sta_average()
    tfnoise_sta = TFNoise(stanoise)
    tfnoise_sta.transfer_func()
    return tfnoise_sta

def test_comply_day_demo():
    daynoise = test_day_average()
    sta = get_meta.get_stdb()
    comply_day = Comply(objnoise=daynoise, sta=sta)
    comply_day.calculate_compliance()
    return comply_day

def test_comply_sta_demo():
    stanoise = test_sta_average()
    sta = get_meta.get_stdb()
    comply_sta = Comply(objnoise=stanoise, sta=sta)
    comply_sta.calculate_compliance()
    return comply_sta
