"""
LDC Solver - Python Prototype
==============================

A 2D incompressible Navier-Stokes solver for lid-driven cavity flow.

This package provides a pure Python implementation using Chorin's 
projection method with explicit diffusion on a uniform staggered grid.

Main modules:
- LDC.py: Main solver implementation
- operators.py: Finite difference operators (gradient, averaging)
- ghia_data.py: Benchmark data from Ghia et al. (1982)
"""

__version__ = "1.0.0"

from . import operators as op
from . import ghia_data as gh

__all__ = ['op', 'gh']
