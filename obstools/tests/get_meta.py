import stdb
import numpy as np
from pkg_resources import resource_filename
from orientpy import utils
from obspy.clients.fdsn import Client


def get_stdb():
    dbfile = resource_filename('orientpy',
                               'examples/data/LOBS3.pkl')
    db = stdb.io.load_db(dbfile)
    return db['YH.LOBS3']


def get_cat():

    sta = get_stdb()
    cat_client = Client()

    tstart = sta.startdate
    tend = sta.enddate
    minmag = 6.0
    try:
        cat = cat_client.get_events(starttime=tstart, endtime=tend,
                                    minmagnitude=minmag, maxdepth=40.)
        reps = np.unique(utils.catclean(cat))
    except:
        raise(Exception("  Fatal Error: Cannot download Catalogue"))

    return cat


