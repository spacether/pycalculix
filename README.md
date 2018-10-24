# pycalculix

pycalculix is a Python 3 library to automate and build finite element analysis (FEA) models in Calculix.

![build passing](https://readthedocs.org/projects/docs/badge/?version=latest)

Website: http://justinablack.com/pycalculix/  
Source Code: https://github.com/spacether/pycalculix  
Documentation: https://pycalculix.readthedocs.io/en/latest/index.html  

#### Usefull applications of Pycalculix
- Trade studies for plane stress, plane strain, or axisymmetric parts
- Quick Kt analysis of 2D geometry
- Learning finite element analyis (FEA) and Python

## Installation
#### Mac OS X
1. Install python3, pycalculix and the fea programs that it uses
```
brew install python3
python3 -mpip install -U numpy
python3 -mpip install -U matplotlib
pip3 install pycalculix
pycalculix-add-feaprograms
```
2. You are done! See 'Usage'

#### Windows
1. Install python3 for [32bit machines](https://www.python.org/ftp/python/3.6.3/python-3.6.3-webinstall.exe) or [64 bit machines](https://www.python.org/ftp/python/3.6.3/python-3.6.3-amd64.exe)
2. In a terminal run the below lines to install needed python libraries, pycalculix, and the fea programs that it uses
```
python -mpip install -U numpy
python -mpip install -U matplotlib
pip install pycalculix
pycalculix-add-feaprograms
```
3. You are done! See 'Usage'

#### Linux, assumes Ubuntu 16.04
1. Install python3-pip, numpy, matplotlib, pycalculix, and the fea programs that it uses
```
sudo apt-get install python3-pip python3-tk
pip3 install --upgrade pip
python3 -mpip install -U numpy
sudo python3 -mpip install -U matplotlib
pip3 install pycalculix
pycalculix-add-feaprograms
```
2. You are done! See 'Sample Program' and 'Usage'

## Sample Program
```
# this is from examples/hole-in-plate-full.py
import pycalculix as pyc

# Vertical hole in plate model, make model
proj_name = 'hole-in-plate-full'
model = pyc.FeaModel(proj_name)
model.set_units('m') # this sets dist units to meters

# Define variables we'll use to draw part geometry
diam = 2.0 # hole diam
ratio = 0.45
width = diam/ratio   #plate width
print('D=%f, H=%f, D/H=%f' % (diam, width, diam/width))
length = 2*width  #plate length

# Draw part geometry, you must draw the part CLOCKWISE, x, y = radial, axial
part = pyc.Part(model)
part.goto(length*0.5, -width*0.5)
part.draw_line_ax(width)
part.draw_line_rad(-length)
part.draw_line_ax(-width)
part.draw_line_rad(length)
hole_lines = part.draw_hole(0, 0, diam*0.5, filled=False)
model.set_ediv(hole_lines, 10)
part.chunk()
model.plot_geometry(proj_name+'_geom') # view the geometry

# set loads and constraints
pressure = -1000
model.set_load('press', part.top, pressure)
model.set_load('press', part.bottom, pressure)
model.set_constr('fix', ['P6', 'P8'], 'y')
model.set_constr('fix', ['P4', 'P7'], 'x')

# set part material
mat = pyc.Material('steel')
mat.set_mech_props(7800, 210*(10**9), 0.3)
model.set_matl(mat, part)

# set the element type and mesh database
model.set_eshape('quad', 2)
model.set_etype('plstress', part, 0.1)
model.mesh(1.0, 'gmsh') # mesh 1.0 fineness, smaller is finer
model.plot_elements(proj_name+'_elem')   # plot part elements
model.plot_pressures(proj_name+'_press')
model.plot_constraints(proj_name+'_constr')

# make and solve the model
prob = pyc.Problem(model, 'struct')
prob.solve()

# view and query results
sx = prob.rfile.get_nmax('Sx')
print('Sx_max: %f' % sx)

# Plot results
fields = 'Sx,Sy,S1,S2,S3,Seqv,ux,uy,utot,ex'    # store the fields to plot
fields = fields.split(',')
for field in fields:
    fname = proj_name+'_'+field
    prob.rfile.nplot(field, fname, display=False)
```

## Usage
1. To run a pycalcuix file you have to have pycalulix installed,
see 'Installation' above
2. The you can then write your own pycalculix programs or run one of the example files here:
https://github.com/spacether/pycalculix/tree/master/examples
3. To run a file:  

#### Windows
- Graphical user interface:
  - Double click the file  
  If the .py extension is associated correctly you can double click it to run the .py program
- Console:
  - cd into the directory with your .py file in it, then in the terminal enter
  ```
  python the_program.py
  ```  
  where the_program.py is the name of the file that you are running
  This assumes that python3 is your active python installation

#### Mac/Linux
- Console:
  - cd into the directory with your .py file in it, then in the terminal enter
  ```
  python3 the_program.py
  ```  
  where the_program.py is the name of the file that you are running


## Elements Supported:
Axisymmetric, plane stress, and plane strain elements are supported.

- First and second order triangles and quadrilaterals are supported.
  - First order elments only have corner nodes
  - Second order elements have corner and midside nodes

Second order elements produce more accurate results
Setting element divisions on lines is supported


## Geometry Building
- One can build separate parts made of points, lines, arcs, and areas.
- Straight lines and arcs are currently supported.
- One can draw a part made of straight lines, then smooth out corners by adding arcs with the part method: part.fillet_lines(L1, L2, arc_radius)
- The new arc will be tangent to both adjacent lines.


## Loading
- Force loading
- Constant pressure
- Linearly varying pressure (water pressure)
- Gravity
- Rotational speed forces
- Displacement constraints
- Loads are stored on geometry primitives (points lines, areas) and
can be applied before or after meshing.


## Files Produced
Meshing and solving are done in the background using cgx or gmsh for
meshing, and Calculix ccx for solving.

Files Used:
- .fbd (Calculix cgx gemetry file)
- .inp (Calculix solver input file, or mesh definition)
- .geo (Gmsh geometry file)
- .msh (Gmsh native mesh file)
- .frd (Calculix ccx nodal results file, values are at nodes and were
  created by interpolating element integration point results back to
  the nodes)
- .dat (Calculix ccx element results file, includes integration point
  results)


## Uninstall
#### Windows
```
pycalculix-remove-feaprograms
pip uninstall pycalculix
```
#### Mac/Linux
```
pycalculix-remove-feaprograms
pip3 uninstall pycalculix
```

## Devlopment
- Download this repo
```
git clone git@github.com:spacether/pycalculix.git
```
- Make python3 virtual env
```
python -m venv venv
```
- Activate it
  - Windows `venv\scripts\activate`
  - Mac/linux `source venv/bin/activate`
- Locally install pycalculix
```
pip install -e .
pycalculix-add-feaprograms
```
- Now any changes that you make to your local version of pycalculix
  will be live in your virtual environment

## License
See LICENSE.txt (Apache 2.0)


## Creator
Justin Black, justin.a.black[at-sign]gmail[dot]com
Initial Release: December 2014


## Change Log

#### 1.0.0
- Fixes the issue where the windows FEA installer hangs:
  https://github.com/spacether/pycalculix/issues/47
- Fixes the issue where the ccx version check in Windows fails
  https://github.com/spacether/pycalculix/issues/48
- Fixes the issue on OS X High Sierra where ccx assumes that the host has gcc7
  installed on it: https://github.com/spacether/pycalculix/issues/50
  - Solution: install gcc@7 and symbolic link it to the expected location


#### 0.9.6
- Fixes part.fillet_lines method, all tests now pass
  - Verified on Mac OS X
  - Closes github Issue 32: '2/15 of the tests fail'

#### 0.9.5 (github + pypi)
- Adds tests: sample tests at tests/test_samples.py
- Adds tests: meshing tests at tests/test_meshing.py
- Adds solving and meshing timeout exception to capture when they hang
- Fixes dxf import feature, syntax updates to use dfxgrabber >= 0.8.0,
  Issue 32
- Adds requirement for dfxgrabber >= 0.8.0 to ensure that dxf import works
- Pegs Mac gmsh version to gmsh == 3.0.5 because version 3.0.6
  throws segault errors when meshing
- Fixes a bug where solver input file does not write material before
  time steps, Issue 32
- Fixed ccx installer on Windows: zip file is now found and
  downloaded
- Throws an exception if ccx version is too old, v 2.7 and earlier is
  too old
- Pegs Win gmsh install version to 3.0.5
- Updates the calculation for element Seqv, S1, S2, and S3 avg max
  and min values. Now calculates Seqv and principal stresses at all
  integration points, then calculates the avg max and min of those
  values
- Win pegged ccx to version 2.12
- Mac brew brewsci/science/calculix-ccx is currently at ccx version 2.13
- Ubuntu apt-get calculix-ccx is currently is currently at version 2.11

#### 0.9.4 (github only)
- removed gmsh and calculix
- moved dist and documentation building and example cleanup into make
  file
- changed the license to Apache 2.0
- added command line tool to install/uninstall gmsh and ccx for
  windows/mac os x/ubuntu
  - pycalculix-add-feaprograms
  - pycalculix-remove-feaprograms
- added requests library requirement for pycalculix-add-feaprograms
- fixed bug where frd files could no longer be read because Calculix
  results keywords changed since initial 2014 release

#### 0.9.3  
- ADDED: multiple parts with contacts
  - See example files: pipe-crush-elastic.py, pinned-plate.py
- ADDED: Import CAD geometry from dxf file
  - See pycalculix.CadImporter
  - Examples:
    - import-dxf.py
- ADDED: Element results plotting added
  - Element results plotting:   pycalculix.Problem.rfile.eplot()
  - Nodal results plotting:     pycalculix.Problem.rfile.nplot()
  - Number Formatting:
    - Strain results now use scientific formatting
    - Others use nearest 10**3 suffixes
  - Max and min values now listed above the colorbar
- ADDED: method to draw an arc by swept angle in degrees
  - part.draw_arc_angle(degrees_ccw, center_x, center_y)
- ADDED: min_val and max_val can now be passed to eplot and nplot
  - This lets the user set the results range that they want to see:
  - min_val <= colored results <= max_val
  - Values under and over are greyed out, darker under, lighter over
- ADDED: internal holes in parts
  - One can make circular holes, or draw complicated holes.
  - See examples:
    - hole-in-plate-full.py, multihole.py
- ADDED: Added set_ediv method to FeaModel class.
  - This method sets the number of elements on a line.
  - line.set_ediv still works.
- ADDED: Robust selection object: feamodel.view
  - This object is feamodel.view Most important methods are:
    - view.select_all, view.select, view.allsel_under
    - All plotting now uses the current selection set
  - SYNTAX: updated how parts, materials, problems are made
    - Make part:
      - pycalculix.FeaModel.make_part or pycalculix.Part
    - Make material:
      - pycalculix.FeaModel.make_matl or pycalculix.Material
    - Make problem (previously called model):
      - pycalculix.FeaModel.make_problem or pycalculix.Problem
    - Make Results File:
      - pycalculix.Problem.rfile or pycalculix.ResultsFile(problem)
- FIX: Plotting fix, closing triangles in the correct direction in
  matplotlib
- DOC: All code separated into modules for clarity.
- DOC: Docstrings added to all classes + methods + functions
- PLOTTING: Closed areas are now filled in yellow when plotting
  geometry.
- PLOTTING: Signed line names are shown and internal to the area.
- BACKEND: Implemented signed line and signed arc class.
  - Pressures can now be applied on these signed lines.
  - Many methods and variables made private to clean up name space.

![Analytics](https://ga-beacon.appspot.com/UA-97855011-1/pycalculix_github?pixel)
