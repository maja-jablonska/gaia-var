from setuptools import setup, find_packages

setup(
    name='gaia_var',
    version='0.9.3',
    packages=find_packages(include=['gaia_var', 'gaia_var.*']),
        install_requires=[
        'astropy~=5.1',
        'pandas~=1.4.4',
        'matplotlib~=3.5.2'
    ],
    description='Gaia variability utility packages',
    long_description='''Package providing utility functions
    such as downloading the lightcurves and variability tables.''',
    author='Maja Jablonska',
    author_email='majajjablonska@gmail.com'
)