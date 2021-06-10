"""MachUpS: A setady-periodic lifting-line simulator."""

from setuptools import setup
import os
import sys

setup(name = 'MachUpS',
    version = '0.0.0',
    description = "MachUpS: A setady-periodic lifting-line simulator.",
    url = 'https://github.com/usuaero/steady-periodic-lifting-line',
    author = 'usuaero',
    author_email = 'doug.hunsaker@usu.edu',
    install_requires = ['numpy>=1.18', 'scipy>=1.4', 'pytest', 'matplotlib'],
    python_requires ='>=3.6.0',
    license = 'MIT',
    packages = ['machupS'],
    zip_safe = False)
