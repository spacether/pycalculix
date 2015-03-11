Description:
pycalculix is a Python 3 library to automate and build finite element analysis (FEA) models in Calculix.
Meshing uses Calculix or GMSH.
Website: http://justinablack.com/pycalculix/
Source Code: https://github.com/spacether/pycalculix
Documentation: http://spacether.github.io/pycalculix/pycalculix.html

Usefull applications of Pycalculix:
-Trade studies for plane stress, plane strain, or axisymmetric parts
-Quick Kt analysis of 2D geometry
-Learning finite element analyis (FEA) and Python

Folder layout:
The folders of this project are laid out to allow me to distribute it on pypi.
Pypi link: https://pypi.python.org/pypi/pycalculix
Note: the main body of pycalculix code is in pycalculix/feamodel.py module

Notes, Cutting Areas:
I built a chunker in python which tries to cut big areas (> 5 sides) which
cgx can't mesh into smaller areas (<= 5 sides) which are meshable in cgx.
The chunker may not always be able to cut areas correctly.

License:
See LICENSE.txt (GPL v2)

Creator:
Justin Black, justin.a.black[at-sign]gmail[dot]com
Initial Release: December 2014

Elements Supported:
Axisymmetric, plane stress, and plane strain elements are supported.
First and second order triangles and quadrilaterals are supported.
  First order elments only have corner nodes
  Second order elements have corner and midside nodes
Second order elements produce more accurate results
Setting element divisions on lines is supported

Geometry Building:
One can build separate parts made of points, lines, arcs, and areas.
Straight lines and arcs are currently supported.
One can draw a part made of straight lines, then smooth out corners by adding
arcs with the part method: part.fillet_lines(L1, L2, arc_radius)
The new arc will be tangent to both adjacent lines.

Loading:
Force loading
Constant pressure
Linearly varying pressure (water pressure)
Gravity 
Rotational speed forces
Displacement constraints
Loads are stored on geometry primitives (points lines, areas) and can be
applied before or after meshing.

Getting Started:
Please read through the below files in the examples folder:
example1_dam.py
example2_hole_in_plate.py
example3_compr_rotor.py
example4_hole_kt.py
example5_times_dam.py

Files Produced:
Meshing and solving are done in the background using cgx or gmsh for meshing,
and Calculix ccx for solving.
Files Used:
*.fbd (Calculix cgx gemetry file)
*.inp (Calculix solver input file, or mesh definition)
*.geo (Gmsh geometry file)
*.msh (Gmsh native mesh file)
*.frd (Calculix ccx nodal results file, values are at nodes and were created
       by interpolating element integration point results back to the nodes)
*.dat (Calculix ccx element results file, includes integration point results)

Installation:
pycalculix requires the below software:
  Python3+
  Numpy (S1,S2,S3 calculation)
  Matplotlib (plotting)
  Calculix (solving)
  Gmsh (meshing)

Easiest Installation:
Install python 3+, then type
pip install pycalculix
Note: This installs the calculix and gmsh programs in sub-folders in the python
      pycalculix folder.
      Running these included binaries from pycalculix only works for windows.
        Linux machines will need libraries to be installed for calculix to work.
	They would also need the programs to be chmoded so they could be run as
	executables.      

Separate Installs:
Install Python 3+: https://www.python.org/downloads/release/python-342/
Install numpy: http://sourceforge.net/projects/numpy/files/NumPy/1.9.1/
Install matplotlib: http://matplotlib.org/downloads.html
Install Calculix: http://www.calculix.de/
    Linux Version: http://www.dhondt.de/
    Windows Version: http://www.bconverged.com/download.php#calculix
Install Gmsh: http://geuz.org/gmsh/#Download
pip install pycalculix
Then pass locations to ccx, cgx, and gmsh per the example on the
pycaculix site: http://justinablack.com/pycalculix/

Anaconda:
	An installation package that includes Python and many python libraries
	Anaconda includes the below Python3+, Numpy, and Matplotlib.
	If you are a Python beginner, I suggest downloading and installing it
	rather than the below separate installers.
	Url:
	http://continuum.io/downloads#py34	

Optional Software:
Suggested IDE (program to edit and run python programs):
Wing IDE:
	http://wingware.com/downloads/wingide-101

Version Updates:
0.9.3
    ADDED: Element results plotting added
        Element results plotting:   pycalculix.Problem.rfile.eplot()
        Nodal results plotting:     pycalculix.Problem.rfile.nplot()
        Number Formatting:
            Strain results now use scientific formatting
            Others use nearest 10**3 suffixes
        Max and min values now listed above the colorbar
    ADDED: method to draw an arc by swept angle in degrees
        part.draw_arc_angle(degrees_ccw, center_x, center_y)
    ADDED: min_val and max_val can now be passed to eplot and nplot
        This lets the user set the results range that they want to see:
        min_val <= colored results <= max_val
        Values under and over are greyed out, darker under, lighter over
    ADDED: internal holes in parts
        One can make circular holes, or draw complicated holes.
        See example 7
    ADDED: Added set_ediv method to FeaModel class.
        This method sets the number of elements on a line.
        line.set_ediv still works.
    ADDED: Robust selection object: feamodel.view
        This object is feamodel.view Most important methods are:
        view.select_all, view.select, view.allsel_under
        All plotting now uses the current selection set
    SYNTAX: updated how parts, materials, problems are made
        Make part: 
            pycalculix.FeaModel.make_part or pycalculix.Part    
        Make material:
            pycalculix.FeaModel.make_matl or pycalculix.Material    
        Make problem (previously called model):
            pycalculix.FeaModel.make_problem or pycalculix.Problem
        Make Results File:
            pycalculix.Problem.rfile or pycalculix.ResultsFile(problem)
    FIX: Plotting fix, closing triangles in the correct direction in matplotlib
    DOC: All code separated into modules for clarity.
    DOC: Docstrings added to all classes + methods + functions
    PLOTTING: Closed areas are now filled in yellow when plotting geometry.
    PLOTTING: Signed line names are shown and internal to the area.
    BACKEND: Implimented signed line and signed arc class.
            Pressures can now be applied on these signed lines.
            Many methods and variables made private to clean up name space.

TODO:
Fix chunking with the real dxf part
    Only 'kontrola.dxf' has problems right now
    Almost there, get rid of zero length curve
Add plotting of internal points
Fix the Linux distribution to make it work out of the box
Add run line to example files for linux
Add tutorial videos
Update pdf
Improve Pylint Scores:
    feamodel:       8.96
    geometry:       8.79
    results_file:   9.75 X
    part:           9.64 X
    selector        9.89 X
    mesh:           9.66 X
    loads:          9.80 X
    cadimporter     9.71 X
    problem:        9.82 X
    base_classes:   9.51 X
    components:     10.0 X
    material:       10.0 X
    environment     10.0 X

Future Goals:
-CAD import of brep and igs via gmsh
-CAD export via gmsh (step, brep)
-Ability to make a new field (% yield etc)
-Double check effective strain calculation:
    http://orange.engr.ucdavis.edu/Documentation12.0/120/ans_thry.pdf 
    pg 23 + 24
    https://support.tnodiana.com/manuals/d944/Analys/node405.html
    https://books.google.com/books?id=70vvzjngyQEC&pg=PA21&lpg=PA21&dq=equivalent+strain+principal+strain&source=bl&ots=bPj4dHfddy&sig=X6MAdyeq34X9uNQRZ1poKXHZK9c&hl=en&sa=X&ei=lr2_VL_lM4q4ggS7loHICw&ved=0CFwQ6AEwCQ#v=onepage&q=equivalent%20strain%20principal%20strain&f=false
-Is my element stress calculation correct? Should it be an average, or a centroid val?
    http://mostreal.sk/html/elem_55/chapter2/ES2-2.htm
    http://orange.engr.ucdavis.edu/Documentation12.0/120/ans_thry.pdf#page=536&zoom=auto,32.4,582.228
-Add area chunker around arcs
-Add part offset function to offset a line and close an area
-Add contact between parts
-Add compression supports
-Reaction force fix?
-Coupling nodes under points and lines ("gluing")
    This would allow 'correct' stresses for knock up at discontinuities
-Add struct-thermal and thermal support
-Order faces under slines to be sequential so they can be output later
-Order nodes under slines to be sequential so they can be output later
-Add mp4 maker
-Bolted joint example perhaps, nodal thickness on bolt and nut areas
-Interactive selection?
-Face plotting?