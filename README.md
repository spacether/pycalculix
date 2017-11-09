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
## Mac OS X
1. Install python3, pycalculix and the fea programs that it uses
```
brew install python3
python3 -mpip install -U numpy
python3 -mpip install -U matplotlib
pip3 install pycalculix
pycalculix-add-feaprograms
```
2. You are done! See 'Usage'

## WINDOWS:
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

## LINUX, assumes Ubuntu 16.04
1. Install python3-pip, numpy, matplotlib, pycalculix, and the fea programs that it uses
```
sudo apt-get install python3-pip
pip3 install --upgrade pip
python3 -mpip install -U numpy
sudo python3 -mpip install -U matplotlib
pip3 install pycalculix
pycalculix-add-feaprograms
```
2. You are done! See 'Usage'


## Usage
1. To run a pycalcuix file you have to have pycalulix installed,
see Installation below
2. The you can then write your own pycalculix programs or run one of the example
files in the eamples folder on github:
https://github.com/spacether/pycalculix/tree/master/examples
3. To run a file:
- WINDOWS:
  - Graphical user interface:
    1. Double click the file: if the .py extension is associated correctly you can double click it to run the .py program
  - Console:
    1. cd into the directory with your .py file in it
    2. type:  
    python the_program.py  
    where the_program.py is the name of the file that you are running
    This assumes that python3 is your active python installation

- LINUX/MAC:
  - Console:
    1. cd into the directory with your .py file in it
    2. type:  
    python3 the_program.py  
    where the_program.py is the name of the file that you are running


## License
See LICENSE.txt (Apache 2.0)


## Creator
Justin Black, justin.a.black[at-sign]gmail[dot]com
Initial Release: December 2014


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




## Change Log

### 0.9.4
- removed gmsh and calculix
- moved dist and documentation building and example cleanup into make file
- changed the license to Apache 2.0
- added command line tool to install gmsh and calculix for: osx

### 0.9.3  
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

TODO:
- Adding the model.view.set_orientation
- convert old view into focus
- Add equation which is coupling
- set_couple('x',['P1', 'P2'])
- set_couple('x',['L1', 'L2'])
- set_couple('x','L1')
- Add tutorial videos
- Update pdf


- Improve Pylint Scores:
  - feamodel        9.28
  - geometry        9.02
  - results_file    9.75
  - part            9.64
  - selector        9.89
  - mesh            9.66
  - loads           9.80
  - cadimporter     9.71
  - problem         9.82
  - base_classes    9.51
  - components      10.0
  - material        10.0
  - environment     9.39
  - connectors      9.57

Future Goals:
- CAD import of: brep, step, iges
- can make geo files with gmsh:
- gmsh freecad_part.iges -o out_iges.geo -0
- gmsh freecad_part.brep -o out_brep.geo -0
- gmsh freecad_part.step -o out_step.geo -0
- Importing the geo file:
- Do points need explicit numbers? no, store them in a dict
- points[geo_str] = point
- Need to convert spline to arcs if it's a circle or arc >= 180
- Will need to delete some points + make center point

- Ability to make a new field (% yield etc)
- Double check effective strain calculation:
    http://orange.engr.ucdavis.edu/Documentation12.0/120/ans_thry.pdf
    pg 23 + 24
    https://support.tnodiana.com/manuals/d944/Analys/node405.html
    https://books.google.com/books?id=70vvzjngyQEC&pg=PA21&lpg=PA21&dq=equivalent+strain+principal+strain&source=bl&ots=bPj4dHfddy&sig=X6MAdyeq34X9uNQRZ1poKXHZK9c&hl=en&sa=X&ei=lr2_VL_lM4q4ggS7loHICw&ved=0CFwQ6AEwCQ#v=onepage&q=equivalent%20strain%20principal%20strain&f=false
- Is my element stress calculation correct? Should it be an average, or a centroid val?
    http://mostreal.sk/html/elem_55/chapter2/ES2-2.htm
    http://orange.engr.ucdavis.edu/Documentation12.0/120/ans_thry.pdf#page=536&zoom=auto,32.4,582.228
- Add area chunker around arcs
- Add part offset function to offset a line and close an area
- Add contact between parts
- Add compression supports
- Reaction force fix?
- Coupling nodes under points and lines ("gluing")
    This would allow 'correct' stresses for knock up at discontinuities
- Add struct-thermal and thermal support
- Order faces under slines to be sequential so they can be output later
- Order nodes under slines to be sequential so they can be output later
- Add mp4 maker
- Bolted joint example perhaps, nodal thickness on bolt and nut areas
- Interactive selection?
- Face plotting?
- Aluminum can model?
  - https://www.youtube.com/watch?v=hUhisi2FBuw

![Justin Analytics](https://ga-beacon.appspot.com/UA-97855011-1/pycalculix_github?pixel)
