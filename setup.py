from setuptools import setup, find_packages

setup(
    name='gaia_toolkit',
    version='0.9.0',
    packages=find_packages(include=['gaia_toolkit', 'gaia_toolkit.*']),
        install_requires=[
        'astropy~=5.1',
        'pandas~=1.4.4',
        'matplotlib~=3.5.2'
    ]
)