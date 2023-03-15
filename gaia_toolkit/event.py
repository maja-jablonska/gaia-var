import matplotlib.pyplot as plt
import pandas as pd
from typing import Dict, List, Optional
from astropy.units import Quantity
from astropy.coordinates import Angle
import astropy.units as u
from gaia_toolkit.archive import (AVAILABLE_COLUMNS, FIVEPARAM_ALL_FIELDS, FIVEPARAM_RA, FIVEPARAM_DEC,
                          fetch_gaia_data, fetch_gaia_photometry, read_fiveparam_fit_from_dataframe)

class GaiaEvent():
    def __init__(self, source_id: int):
        self.__source_id = source_id
        self.__5param_fit: Dict[str, Optional[float]] = {f: None for f in FIVEPARAM_ALL_FIELDS}
        self.__lightcurve: Optional[pd.DataFrame] = None
        self.__gaia_archive_data: Optional[pd.DataFrame] = None

    @property
    def source_id(self) -> int:
        return self.__source_id

    @property
    def fiveparam_fit(self) -> Dict[str, Optional[Quantity]]:
        return self.__5param_fit

    @property
    def photometry(self) -> pd.DataFrame:
        return self.__lightcurve

    @property
    def ra(self) -> Optional[float]:
        return self.__5param_fit[FIVEPARAM_RA]

    @property
    def dec(self) -> Optional[float]:
        return self.__5param_fit[FIVEPARAM_DEC]
    
    @property
    def coordinates_string(self) -> str:
        ra = self.ra
        dec = self.dec
        if ra and dec:
            ra_angle: Angle = Angle(ra, unit=u.deg)
            dec_angle: Angle = Angle(dec, unit=u.deg)
            return f'{ra_angle.to_string(sep=":")} {dec_angle.to_string(sep=":")}'
        else:
            return ""

    @property
    def reference_epoch(self) -> float:
        return self.__reference_epoch

    @property
    def fiveparam_fit_values(self) -> Dict[str, Optional[float]]:
        return {f: k.value for f, k in self.__5param_fit.items() if k is not None}

    @property
    def gaia_archive_data(self) -> Dict[str, Optional[pd.DataFrame]]:
        return self.__gaia_archive_data

    def band_photometry(self, band: str) -> pd.DataFrame:
        return self.photometry[self.photometry['band']==band]

    @photometry.setter
    def photometry(self, new_photometry: pd.DataFrame):
        self.__lightcurve = new_photometry

    @reference_epoch.setter
    def reference_epoch(self, new_reference_epoch: float):
        self.__reference_epoch = new_reference_epoch

    @fiveparam_fit.setter
    def fiveparam_fit(self, new_fiveparam: Dict[str, Optional[float]]):
        self.__5param_fit = new_fiveparam

    def fetch_gaia_archive_data(self):
        self.__gaia_archive_data: pd.DataFrame = fetch_gaia_data(self.source_id, AVAILABLE_COLUMNS)
        self.fiveparam_fit = read_fiveparam_fit_from_dataframe(self.__gaia_archive_data)

    def fetch_dr3_photometry(self):
        self.photometry = fetch_gaia_photometry(self.source_id)


    def plot_photometry(self, ax: Optional[plt.Axes] = None, bands: Optional[List[str]] = None):
        BAND_COLORS: Dict = {'G': 'slategray', 'BP': 'royalblue', 'RP': 'firebrick'}
        
        if not ax:
            _, ax = plt.subplots(1, 1)

        unique_bands = self.photometry['band'].unique()

        if not bands:
            bands = unique_bands
        elif not all([b in unique_bands for b in bands]):
            raise ValueError(f'Requested bands not present in data. Select bands from {", ".join(unique_bands)}')

        for band in bands:
            band_photometry = self.band_photometry(band)
            ax.errorbar(band_photometry.time_year,
                        band_photometry.mag,
                        yerr=band_photometry.mag_err,
                        fmt='o', ms=2, elinewidth=.5, color=BAND_COLORS.get(band, 'black'), label=band)

        ax.set_xlabel('time [year]')
        ax.set_ylabel('luminosity [mag]')
        ax.invert_yaxis()
        ax.legend()
