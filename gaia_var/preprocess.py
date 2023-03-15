from astropy.time import Time
import numpy as np
import pandas as pd


def gaia_time_to_bjd(tcb_time: float) -> Time:
    """Convert Gaia onboard time to JD time

    Args:
        tcb_time (float): Gaia onboard time

    Returns:
        Time: astropy Time object
    """
    return Time(tcb_time + 2455197.5, format='jd')
