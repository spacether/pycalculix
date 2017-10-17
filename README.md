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


## Folder layout
The folders of this project are laid out to allow me to distribute it on pypi.
Pypi link: https://pypi.python.org/pypi/pycalculix

Note: the main body of pycalculix code is in pycalculix/feamodel.py module


## Notes, Cutting Areas
I built a chunker in python which tries to cut big areas (> 5 sides) which
cgx can't mesh into smaller areas (<= 5 sides) which are meshable in cgx.
The chunker may not always be able to cut areas correctly.


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


## Getting Started
1. To run a pycalcuix file you have to have pycalulix installed,
see Installation below
2. The you can then write your own pycalculix programs or run one of the example
files in the eamples folder on github:
https://github.com/spacether/pycalculix/tree/master/examples
3. To run a file:
  ##### WINDOWS:

  - Graphical user interface:
    1. Double click the file: if the .py extension is associated correctly you can double click it to run the .py program

  - Console:
    1. cd into the directory with your .py file in it
    2. type:  
    python the_program.py  
    where the_program.py is the name of the file that you are running
    This assumes that python3 is your active python installation

  ##### LINUX:

  - Console:
    1. cd into the directory with your .py file in it
    2. type:  
    python3 the_program.py  
    where the_program.py is the name of the file that you are running


## Files Produced
Meshing and solving are done in the background using cgx or gmsh for meshing, and Calculix ccx for solving.

Files Used:
- .fbd (Calculix cgx gemetry file)
- .inp (Calculix solver input file, or mesh definition)
- .geo (Gmsh geometry file)
- .msh (Gmsh native mesh file)
- .frd (Calculix ccx nodal results file, values are at nodes and were created by interpolating element integration point results back to the nodes)
- .dat (Calculix ccx element results file, includes integration point results)


## Installation
pycalculix requires the below software:
- Python3+
- Numpy (S1,S2,S3 calculation)
- Matplotlib (plotting)
- Calculix (solving)
- Gmsh (meshing)

## Mac OS X
1. Install python3, pycalculix and the fea programs that it uses
```
brew install python3
pip3 install pycalculix
pycalculix-add-feaprograms
```
2. You are done! See 'Getting Started'

##### WINDOWS:
1. Install Anaconda* python 3.4:  
  http://continuum.io/downloads#py34  
  This includes required libraries like numpy and matplotlib.

  Note (32 vs 64 bit):

  - 64-bit:

    This is an easier option for 64 bit systems because the binaries of the other libraries are harder to find.

  - 32-bit:

    Rather than installing anaconda you can just install python 3+ from the python site:  

    https://www.python.org/downloads/release/python-343/

    When you install pycalculix in the next step below, the required libraries will autoinstall on your system.

2. In a console window, type:  
```
pip install pycalculix
pycalculix-add-feaprograms
```
Note: the second line installs the calculix and gmsh programs on your computer
Running these included binaries from pycalculix only works for windows.

3. You are done! See 'Getting Started'

##### LINUX
(assumes Ubuntu 14.04)

1. Install required prerequsisites. In the console enter:  
  sudo apt-get install python3-pip python3-matplotlib gmsh  
  Note: this installs
  - pip (a python library downloader) for python 3  
  - python3 matplotlib (the required matplotlib library needed for plotting)  
  - gmsh (software needed to mesh you FEA models)

2. Install both calculix ccx and calculix cgx for your architecture (32 or 64 bit)

  Ubuntu 14.04.1 (Trusty)  
  32-bit:  
  ccx:   https://code.launchpad.net/~cae-team/+archive/ubuntu/ppa/+build/7043228/+files/calculix-ccx_2.7-0%7E1%2B6%7Eubuntu14.04.1_i386.deb  
  cgx:   https://code.launchpad.net/~cae-team/+archive/ubuntu/ppa/+build/7043230/+files/calculix-cgx_2.7-0%7E1%2B3%7Eubuntu14.10.1_i386.deb

  64-bit:  
  ccx:   https://code.launchpad.net/~cae-team/+archive/ubuntu/ppa/+build/7043227/+files/calculix-ccx_2.7-0%7E1%2B6%7Eubuntu14.04.1_amd64.deb  
  cgx:   https://code.launchpad.net/~cae-team/+archive/ubuntu/ppa/+build/7043229/+files/calculix-cgx_2.7-0%7E1%2B3%7Eubuntu14.10.1_amd64.deb

3. In the console window type:  
pip3 install pycalculix

4. You are done! See 'Getting Started'


## Separate Installs
Install Python 3+: https://www.python.org/downloads/release/python-342/  
Install numpy: http://sourceforge.net/projects/numpy/files/NumPy/1.9.1/  
Install matplotlib: http://matplotlib.org/downloads.html  
Install Calculix: http://www.calculix.de/  
- Linux Version: http://www.dhondt.de/
- Windows Version: http://www.bconverged.com/download.php#calculix
Install Gmsh: http://geuz.org/gmsh/#Download  
Install pycalculix: pip install pycalculix  
Then pass locations to ccx, cgx, and gmsh per the example on the  
pycaculix site: http://justinablack.com/pycalculix/  

## Anaconda
An installation package that includes Python and many python libraries
Anaconda includes the below Python3+, Numpy, and Matplotlib.
If you are a Python beginner, I suggest downloading and installing it
rather than the separate installers.
http://continuum.io/downloads#py34

Optional Software:
Suggested IDE (program to edit and run python programs):
Wing IDE:
	http://wingware.com/downloads/wingide-101

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
