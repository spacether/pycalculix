# pycalculix

pycalculix is a Python 3 library to automate and build finite element analysis (FEA) models in Calculix.

Meshing uses Calculix or GMSH
- Website: http://justinablack.com/pycalculix/  
- Source Code: https://github.com/spacether/pycalculix  
- Documentation: http://pythonhosted.org/pycalculix/

Usefull applications of Pycalculix:
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
Note: the second line installs the calculix and gmsh programs on your computer
Running these included binaries from pycalculix only works for windows.  
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
2. You are done! See 'Usage'


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
- Loads are stored on geometry primitives (points lines, areas) and can be
applied before or after meshing.


## Files Produced
Meshing and solving are done in the background using cgx or gmsh for meshing, and Calculix ccx for solving.

Files Used:
- .fbd (Calculix cgx gemetry file)
- .inp (Calculix solver input file, or mesh definition)
- .geo (Gmsh geometry file)
- .msh (Gmsh native mesh file)
- .frd (Calculix ccx nodal results file, values are at nodes and were created by interpolating element integration point results back to the nodes)
- .dat (Calculix ccx element results file, includes integration point results)


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


## License
See LICENSE.txt (Apache 2.0)


## Creator
Justin Black, justin.a.black[at-sign]gmail[dot]com
Initial Release: December 2014


## Change Log

#### 0.9.4 (github only)
- removed gmsh and calculix
- moved dist and documentation building and example cleanup into make file
- changed the license to Apache 2.0
- added command line tool to install/uninstall gmsh and ccx for windows/mac os x/ubuntu
  - pycalculix-add-feaprograms
  - pycalculix-remove-feaprograms

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
- FIX: Plotting fix, closing triangles in the correct direction in matplotlib
- DOC: All code separated into modules for clarity.
- DOC: Docstrings added to all classes + methods + functions
- PLOTTING: Closed areas are now filled in yellow when plotting geometry.
- PLOTTING: Signed line names are shown and internal to the area.
- BACKEND: Implemented signed line and signed arc class.
  - Pressures can now be applied on these signed lines.
  - Many methods and variables made private to clean up name space.

## Developer Todo:
- fix the error where FRD files can no longer be read
- confirm that set_ediv still works, remove the past merge if that broke it
