from collections import defaultdict
from astroquery.gaia import Gaia
from typing import Dict, List
import pandas as pd
import astropy.units as u
import numpy as np
import pandas as pd
from typing import List, Optional
from ssl import SSLEOFError

from gaia_var.preprocess import gaia_time_to_bjd

GAIA_TABLE_DR3 = pd.DataFrame()
AVAILABLE_COLUMNS = []

try:
    GAIA_TABLE_DR3 = Gaia.load_table('gaiadr3.gaia_source')
    AVAILABLE_COLUMNS = [t.name for t in GAIA_TABLE_DR3.columns]
except SSLEOFError:
    print("Couldn't connect to Gaia server! check the internet connection and try again...")
except Exception as e:
    print(f'Unknown exception occured while connecting to Gaia Archive: {e}')

MAX_RANDOM_INDEX: int = 1811709770


COLUMNS_TO_FETCH: List[str] = [
    'solution_id', 'source_id',
    'ra', 'ra_error', 'dec', 'dec_error',
    'parallax', 'parallax_error',
    'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
    'ruwe', 'phot_g_mean_mag', 'phot_rp_mean_mag', 'phot_bp_mean_mag'
    ]


CORRELATION_COLUMNS: List[str] = [
    'ra_dec_corr',
    'ra_pmra_corr', 'ra_pmdec_corr', 'dec_pmra_corr', 'dec_pmdec_corr',
    'ra_parallax_corr', 'dec_parallax_corr', 'parallax_pmra_corr', 'parallax_pmdec_corr',
    'pmra_pmdec_corr'
    ]


METADATA_COLUMNS: List[str] = [
    'phot_variable_flag', 'has_xp_sampled', 'has_rvs',
    'has_epoch_photometry', 'non_single_star', 'has_epoch_rv',
    'in_galaxy_candidates', 'in_qso_candidates', 'in_andromeda_survey'
    ]


QUALITY_COLUMNS: List[str] = [
    'astrometric_n_good_obs_al', 'astrometric_excess_noise', 'astrometric_chi2_al',
    'pseudocolour', 'pseudocolour_error', 'duplicated_source',
    'phot_rp_n_blended_transits', 'phot_bp_n_blended_transits'
    ]

COLOR_COLUMNS: List[str] = [
    'bp_rp'
    ]

ALL_ADDITIONAL_COLUMNS = CORRELATION_COLUMNS + METADATA_COLUMNS + QUALITY_COLUMNS + COLOR_COLUMNS


FIVEPARAM_RA: str = 'ra'
FIVEPARAM_DEC: str = 'dec'
FIVEPARAM_PMRA: str = 'pmra'
FIVEPARAM_PMRA_ERR: str = 'pmra_error'
FIVEPARAM_PMDEC: str = 'pmdec'
FIVEPARAM_PMDEC_ERR: str = 'pmdec_error'
FIVEPARAM_PARALLAX: str = 'parallax'
FIVEPARAM_PARALLAX_ERR: str = 'parallax_error'

FIVEPARAM_PARAMETERS: List[str] = [FIVEPARAM_RA, FIVEPARAM_DEC, FIVEPARAM_PMRA, FIVEPARAM_PMDEC, FIVEPARAM_PARALLAX]
DEFAULT_FIVEPARAM_PARAMETER_UNITS: Dict = {
    FIVEPARAM_RA: u.deg,
    FIVEPARAM_DEC: u.deg,
    FIVEPARAM_PMRA: u.mas/u.year,
    FIVEPARAM_PMDEC: u.mas/u.year,
    FIVEPARAM_PARALLAX: u.mas}

FIVEPARAM_ERRORS: List[str] = [FIVEPARAM_PMRA_ERR, FIVEPARAM_PMDEC_ERR, FIVEPARAM_PARALLAX_ERR]

FIVEPARAM_ALL_FIELDS: List[str] = FIVEPARAM_PARAMETERS + FIVEPARAM_ERRORS

def read_fiveparam_fit_from_dataframe(df: pd.DataFrame, units: Dict = DEFAULT_FIVEPARAM_PARAMETER_UNITS) -> Dict[str, float]:
    if not all([(f in df.columns) for f in FIVEPARAM_ALL_FIELDS]):
        raise ValueError(f'Dataframe doesn\'t contain all fiveparameter fit parameters: {", ".join(FIVEPARAM_ALL_FIELDS)}')

    if len(df.index)<1:
        raise ValueError(f'Dataframe can\'t be empty!')
    
    return {f: (df[f].values[0]*units[f.replace('_error', '')]).to(DEFAULT_FIVEPARAM_PARAMETER_UNITS[f.replace('_error', '')]) for f in FIVEPARAM_ALL_FIELDS}


def mag_error(flux_over_error: float) -> float:
    return 1/(flux_over_error*2.5/np.log(10))


def _fetch_gaia_datalink_product(source_id: int, product_type: str) -> pd.DataFrame:
    retrieval_type = product_type
    data_structure = 'INDIVIDUAL'
    data_release   = 'Gaia DR3'
    datalink = Gaia.load_data(
        ids=[source_id],
        data_release=data_release,
        retrieval_type=retrieval_type,
        data_structure=data_structure,
        verbose=True,
        output_file=None)
    dl_keys  = [inp for inp in datalink.keys()]

    if len(dl_keys)>0:
        print(f'The following Datalink product has been downloaded:')
        for dl_key in dl_keys:
            print(f' * {dl_key}')
    else:
        print(f'No Datalink product downloaded.')
        return pd.DataFrame()
    
    try:
        table: pd.DataFrame = datalink[dl_keys[0]][0].to_table().to_pandas()
    except Exception as e:
        print(f'Error: {e}')
        return pd.DataFrame
    
    return table


def fetch_gaia_photometry(source_id: int) -> pd.DataFrame:
    table =  _fetch_gaia_datalink_product(source_id, 'EPOCH_PHOTOMETRY')
    try:
        table['time_year'] = gaia_time_to_bjd(table['time']).jyear
        table['mag_err'] = mag_error(table['flux_over_error'])
    except Exception as e:
        print(f'Error: {e}')
        return table
    return table


def fetch_gaia_rvs(source_id: int) -> pd.DataFrame:
    return _fetch_gaia_datalink_product(source_id, 'RVS')


def fetch_gaia_xp(source_id: int) -> pd.DataFrame:
    return _fetch_gaia_datalink_product(source_id, 'XP_SAMPLED')


def __source_id_condition(source_id: np.ndarray) -> str:
    source_id_condition: str = ''
    if isinstance(source_id, int):
        print("Fetching for single source...")
        source_id_condition = '= {:d}'.format(source_id)
    elif isinstance(source_id, list) or isinstance(source_id, np.ndarray):
        print(f'Fetching for {len(source_id)} sources...')
        source_id_condition = f'IN ({",".join(["{:d}".format(int(si)) for si in source_id])})'
    else:
        raise ValueError('Source ID must be a single integer, a list, or ndarray of integers!')
    
    return source_id_condition


def fetch_gaia_data(source_id: np.ndarray, additional_columns: Optional[List[str]] = None) -> pd.DataFrame:
    
    source_id_condition = __source_id_condition(source_id)

    additional_columns = additional_columns if additional_columns else []

    column_availability = [(ac in AVAILABLE_COLUMNS) for ac in additional_columns]
    if not all(column_availability):
        additional_columns = np.array(additional_columns)
        column_availability = np.array(column_availability)
        raise ValueError('Columns ' + ', '.join([str(a) for a in additional_columns[~column_availability]]) + ' not in gaiadr3.gaia_source!')

    columns = COLUMNS_TO_FETCH + additional_columns

    try:
        job = Gaia.launch_job_async(f"""select {','.join(set(columns))}
        from gaiadr3.gaia_source where source_id {source_id_condition})
        """)

        results = job.get_results().to_pandas()

        print(f'Fetched {len(results.index)} sources.')

        return results
    except Exception as e:
        print('Exception while querying sources' + getattr(e, 'message', repr(e)))

        return pd.DataFrame(data=[], columns=AVAILABLE_COLUMNS)


VARIABILITY_CLASSES = [
    'rrlyrae',
    'cepheid',
    'planetary_transit',
    'short_timescale',
    'long_period_variable',
    'eclipsing_binary',
    'rotation_modulation',
    'ms_oscillator',
    'agn',
    'microlensing',
    'compact_companion'
]


def get_variability_class(source_id: np.ndarray) -> Optional[str]:

    source_id_condition = __source_id_condition(source_id)

    job = Gaia.launch_job_async(f"""select source_id, phot_variable_flag
    from gaiadr3.gaia_source where source_id {source_id_condition})
    """)

    try:
        results = job.get_results().to_pandas()
        print(f'Fetched variability flags for {len(results.index)} sources.')
    except Exception as e:
        print('Exception while querying phot_variability_flag for sources' + getattr(e, 'message', repr(e)))
        return pd.DataFrame(data=[], columns=AVAILABLE_COLUMNS)

    variable_sources = results[results['phot_variable_flag']=='VARIABLE']
    variable_source_ids = variable_sources['source_id'].values

    variable_source_id_condition = __source_id_condition(variable_source_ids)

    try:
        print(f"""select source_id, {','.join([f'in_vari_{vc}' for vc in VARIABILITY_CLASSES])}
        from gaiadr3.vari_summary where source_id {variable_source_id_condition})
        """)
        job = Gaia.launch_job_async(f"""select source_id, {','.join([f'in_vari_{vc}' for vc in VARIABILITY_CLASSES])}
        from gaiadr3.vari_summary where source_id {variable_source_id_condition})
        """)

        results = job.get_results().to_pandas()

        print(f'Fetching for variability classes for {variable_source_ids}...')

    except Exception as e:
        print('Exception while querying for variability classes for sources' + getattr(e, 'message', repr(e)))
        return pd.DataFrame(data=[], columns=AVAILABLE_COLUMNS)
        
    all_results = defaultdict(lambda: defaultdict(list))

    for _, row in results.iterrows():
        for vc in VARIABILITY_CLASSES:
            if row[f'in_vari_{vc}']:
                job = Gaia.launch_job_async(f"""select * 
                from gaiadr3.vari_{vc} where source_id={row['source_id']})
                """)
                source_result = job.get_results().to_pandas()
                all_results[row['source_id']][vc].append(source_result)

    return all_results
