import pandas as pd
from astropy.time import Time
from gaia_var.preprocess import gaia_time_to_bjd
import os

gaps = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/astrometry_EDR3_gaps.csv'))
gaps['start_jyear'] = gaia_time_to_bjd(gaps['start']).jyear
gaps['end_jyear'] = gaia_time_to_bjd(gaps['end']).jyear

decontamination = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/decontamination.csv'))

def gap_for_time(time: Time) -> pd.DataFrame:
    """Return the gap(s_ for given time

    Args:
        time (Time): astropy Time in question

    Returns:
        pd.DataFrame: gaps occuring during the passed time
    """
    return gaps[(gaps.start_jyear<time.jyear) & (gaps.end_jyear>time.jyear)]

def is_in_gap(time: Time) -> bool:
    """
    Args:
        time (Time): astropy Time in question

    Returns:
        bool: did a gap occur during the passed time
    """
    return len(gap_for_time(time).index)>0

def get_decontamination_times() -> pd.DataFrame:
    return decontamination[decontamination.event!="Refocus"]

def get_refocus_times() -> pd.DataFrame:
    return decontamination[decontamination.event=="Refocus"]
