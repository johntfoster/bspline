# -*- coding: utf-8 -*-
#
"""Compute B-spline basis functions via Cox - de Boor algorithm.

Submodules:
    bspline.bspline
        OO interface (class Bspline)
    bspline.splinelab
        MATLAB-style interface and helper functions.

By default, the Bspline class from bspline.bspline is imported into this namespace when this module is loaded.
"""

from __future__ import absolute_import

# This is extracted automatically by the top-level setup.py.
__version__ = '0.1.1'

# add any imports here, if you wish to bring things into the library's top-level namespace when the library is imported.
from .bspline import Bspline

