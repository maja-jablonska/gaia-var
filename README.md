# Gaia Variables Toolkit
Gaia utilities toolkit for variable events.

## Tools

### Setup

Only Gaia DR3 source ID is required.

```python
import gaia_toolkit as gt
SOURCE_ID: int = 6059400613544951552
event: GaiaEvent = gt.GaiaEvent(SOURCE_ID)
```

### Fetching lightcurves

```python
event.fetch_dr3_photometry()
```

Photometry can be then accessed as a Pandas DataFrame:

```python
event.photometry
```

as well as individual band:

```python
event.band_photometry("G")
```

and can be automatically plotted:
```python
# bands is an optional argument - when it's not passed, all passbands are plotted
event.plot_photometry(bands=["G"])
```

### 5-parameter fit

Can be accessed in a similar way to photometry, returning a dictionary of astroPy Quantities:

```python
event.fetch_gaia_archive_data()
event.fiveparam_fit
```

as well as individual values, but always make note of the units used:

```python
event.fiveparam_fit_values
```

## TO-DOs

### BP/RP spectra

### Variability tables
