from setuptools import setup, find_packages
import os

progs = ['gmsh', 'calculix']
oses = ['win32', 'win64', 'linux32', 'linux64']
folders = []
for program in progs:
    for o in oses:
        folders.append(program+'_'+o)

# platform specific builds
# http://stackoverflow.com/questions/6469508/is-it-possible-to-express-a-platform-specific-dependency-in-setup-py-without-bui
# check -format
# check bdst -format=
# check -p for platform
# https://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/
# version suggestions
# http://semver.org/

setup(
    name = "pycalculix",
    install_requires = ['matplotlib >= 1.3.1', 'numpy'],
    version = "0.9.3",
    description = "Python 3 library to automate and build finite element analysis (FEA) models in Calculix.",
    author = "Justin Black",
    author_email = "justin.a.black@gmail.com",
    packages = find_packages(),
    include_package_data=True,
    url = "http://justinablack.com/pycalculix/",
    keywords = ["FEA", "Finite Element Analysis", "Calculix", "Mechanical Engineering", "CAD"],
    classifiers = [
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering"],
    long_description = """\
'''
Python 3 library to automate and build finite element analysis (FEA) models in Calculix.
----------------------------------------------------------------------------------------
Meshing uses Calculix or GMSH.
Website: http://justinablack.com/pycalculix/
Source Code: https://github.com/spacether/pycalculix

Usefull applications of Pycalculix:
-Trade studies for plane stress, splane strain, or axisymmetric parts
-Quick Kt analysis of 2D geometry
-Leaning finite element analyis (FEA) and Python

Notes:
I build a chunker in python which tries to cut big areas (> 5 sides) which
cgx can't mesh into smaller areas (<= 5 sides) which are meshable in cgx.
The chunker may not always be able to cut areas correctly.

Elements Supported:
Axisymmetric, plane stress, and plane strain elements are supported.
First and second order triangles and quadrilaterals are supported.
  First order elments only have corner nodes
  Second order elements have midside nodes
Second order elements produce more accurate results
Setting element divisions on lines is supported

Geometry Building:
One can build separate parts made of points, lines, arcs, and areas.
Straight lines and arcs are currently supported.
One can draw a part made of straight lines, then smooth out corners by adding
blends/fillets with the part method: part.fillet_lines(L1, L2, arc_radius)
The new filleet will be tangent to both adjacent lines.

Loading:
Force loading, constant pressure, linearly varying pressure, gravity, and rotational speed forces are implemented.
Displacement constraints are also supported.
Loads are stored on geometry primitives (points lines etc) and can be applied before or after meshing.

Examples:
Visit: http://justinablack.com/pycalculix/ to browse example models.
"""
)