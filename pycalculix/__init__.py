#!/usr/bin/python
"""
pycalculix Description:
This module provides a python interface to the calculix FEA preprocessor and
solver. Simplicity is emphasized.

2D models are supported including:
plane stress, plane strain, and axisymmetric problems

Library Capabilities:
Build Geometry
Mesh Geometry
Apply Loads
Solve Model
View Results
Query Results
Select Sub-Portion of Model and view those results

See documentation and examples at:
http://justinablack.com/pycalculix/
https://github.com/spacether/pycalculix
"""

__author__ = "Justin Black"
__copyright__ = "Copyright 2014, Justin Black"
__credits__ = ["Justin Black"]
__license__ = "GPL V2"
__version__ = "0.9.3"
__maintainer__ = "Justin Black"
__email__ = "justin.a.black[insert~at~sine]gmail.com"
__status__ = "Beta"
__all__ = ['environment', 'geometry', 'components', 'loads', 'part', 'mesh',
           'material', 'feamodel', 'model', 'results_file']

from .feamodel import FeaModel
