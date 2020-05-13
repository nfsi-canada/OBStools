from obstools.atacr import DayNoise, StaNoise, TFNoise, EventStream

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

def test_tfnoise_sta_demo():
    stanoise = test_sta_average()
    tfnoise_day = TFNoise(stanoise)
    tfnoise_day.transfer_func()

