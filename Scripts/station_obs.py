import pickle
from obspy.core import UTCDateTime

# Dictionary with all station information
# Keys are station names, and each item contains 
# a list of network names, station names, station lat,
# station lon, and channel type (either BH or HH)
stations = \
        {'FN02C':['7D', 'FN02C', 46.9497, -124.428, 67., 141., 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-28')],\
        'FN03C':['7D', 'FN03C', 46.8872, -124.5251, 93., 26., 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-28')],\
        'FN04C':['7D', 'FN04C', 46.9175, -124.6013, 104., 163., 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-28')],\
        'FN05C':['7D', 'FN05C', 46.8575, -124.6556, 123., 357., 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-25')],\
        'FN07C':['7D', 'FN07C', 46.8554, -124.7860, 158., 171, 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-28')],\
        'FN07A':['7D', 'FN07A', 46.8555, -124.7865, 154., 118, 'HH', UTCDateTime('2011-07-28'), UTCDateTime('2012-07-20')],\
        'FN08C':['7D', 'FN08C', 46.8554, -124.8760, 176., 324, 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-28')],\
        'FN08A':['7D', 'FN08A', 46.8888, -124.8769, 177., 128, 'HH', UTCDateTime('2011-07-29'), UTCDateTime('2012-07-22')],\
        'FN11C':['7D', 'FN11C', 46.8205, -125.0454, 619., 252, 'HH', UTCDateTime('2013-08-31'), UTCDateTime('2014-06-28')],\
        'FN12A':['7D', 'FN12A', 46.8885, -125.1192, 650., 312, 'HH', UTCDateTime('2011-07-27'), UTCDateTime('2012-07-20')],\
        'FN12C':['7D', 'FN12C', 46.8887, -125.1190, 656., 318, 'HH', UTCDateTime('2013-09-03'), UTCDateTime('2014-06-30')],\
        'FN14A':['7D', 'FN14A', 46.0248, -124.9647, 173., 347, 'HH', UTCDateTime('2011-07-31'), UTCDateTime('2012-07-20')],\
        'FN14C':['7D', 'FN14C', 47.0249, -124.9642, 173., 221, 'HH', UTCDateTime('2011-07-31'), UTCDateTime('2012-07-20')],\
        'J52C':['7D', 'J52C', 46.992, -127.0158, 2640, 129., 'BH', UTCDateTime('2013-08-20'), UTCDateTime('2014-05-31')],\
        'J30A':['7D', 'J30A', 45.4242, -128.9069, 2824, 69, 'BH', UTCDateTime('2011-11-20'), UTCDateTime('2012-05-20')],\
        'J30C':['7D', 'J30C', 45.4261, -128.9100, 2786, 65, 'BH', UTCDateTime('2013-08-07'), UTCDateTime('2014-05-15')],\
        'J28A':['7D', 'J28A', 45.0636, -127.1564, 2867, 206, 'BH', UTCDateTime('2011-11-16'), UTCDateTime('2012-05-16')],\
        'J28B':['7D', 'J28B', 45.0631, -127.1552, 2866, 23, 'BH', UTCDateTime('2012-08-24'), UTCDateTime('2013-06-04')],\
        'J28C':['7D', 'J28C', 45.0617, -127.1567, 2889, 186, 'BH', UTCDateTime('2013-08-02'), UTCDateTime('2014-05-21')],\
        'J39C':['7D', 'J39C', 46.18, -129.64, 2656, 268, 'BH', UTCDateTime('2013-08-06'), UTCDateTime('2014-05-16')],\
        }

pickle.dump(stations, open('stations_obs.pkl', 'wb'))

