Description:
pycalculix is a Python 3 library to automate and build finite element analysis (FEA) models in Calculix.
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

License:
See LICENSE.txt (GPL v2)

Creator:
Justin Black, justin.a.black[at-sine]gmail[dot]com
Initial Release: December 2014

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

Getting Started:
Please read through the documented example files below:
example1_dam.py
example2_hole_in_plate.py
example3_compr_rotor.py
example4_hole_kt.py
example5_times_dam.py

Files Produced:
Meshing and solving are done in the background using cgx or gmsh for meshing, and Calculix ccx for solving.
Files Used:
*.fbd (Calculix cgx gemetry file)
*.inp (Calculix solver input file, or mesh definition)
*.geo (Gmsh geometry file)
*.msh (Gmsh native mesh file)
*.frd (Calculix ccx nodal results file, values are at nodes and were created by interpolating element integration point results back to the nodes)
*.dat (Calculix ccx element results file, includes integration point results)

Installation:
Install the below required software.
Then move pycalculix.py to the directory you want to work in.
Any new models should be created in separate .py files. See the above example files in the 'Getting Started' section.

Required Software:
Calculix:
	Free and very full-featured preprocessor (cgx) and finite element analysis (FEA) solver (ccx)
	http://www.calculix.de/
		Linux Version: http://www.dhondt.de/
		Windows Version: http://www.bconverged.com/download.php#calculix
Gmsh: 
	Free and full featured mesher
	http://geuz.org/gmsh/#Download

Anaconda:
	An installation package that includes Python and many often used python libraries
	Anaconda includes the below Python3+, Numpy, and Matplotlib.
	If you are a Python beginner, I suggest downloading and installing it rather than the below separate installers.
	Url:
	http://continuum.io/downloads#py34
	
Python 3+:
	https://www.python.org/downloads/release/python-342/
Numpy: 
	http://sourceforge.net/projects/numpy/files/NumPy/1.9.1/
Matplotlib: 
	http://matplotlib.org/downloads.html

Optional Software:
Suggested IDE (program to edit and run python programs):
Wing IDE:
	http://wingware.com/downloads/wingide-101
	
Future Goals:
-Plotting
  Add element results plotting
-Add compression supports
-Make lists for lines and signed lines (need to write a signed lines class)
-Add struct-thermal and thermal support
-Auto-detect contact regions between parts
-CAD import of brep and igs via gmsh
-CAD export via gmsh (step, brep)
-Add mp4 maker
-Bolted joint example perhaps, nodal thickness on bolt and nut areas
-Interactive selection?