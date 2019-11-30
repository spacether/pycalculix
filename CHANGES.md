## Change Log

#### 1.1.3
- Removes gmsh exerimental flag Mesh.CharacteristicLengthFromCurvature
- Fixes typo in geometry.py

#### 1.1.2
- Sets part.left/right/top/bottom using geometry.ACC constant, tests added
- Fixes issue https://github.com/spacether/pycalculix/issues/57 where part.left
  was not being set if a line was slightly skewed
- Fixes a mac install issue where gcc@7 was assuming a specific X.X.X version
  of gcc7

#### 1.1.1
- Omits test_pinned_plate from Mac OS X with Python >= 3.6 because it does not
converge in ccx

#### 1.1.0
- Adds Microsoft Azure continuous integration tests to verify that the library is working. The following environments are tested:
  - Windows Server 2016 (x64 and x86 architecture)
  - Mac OS X (x64)
  - Ubuntu 16.04 (x64)
- Fixes a bug where cad files could not be loaded if they had a period in their path

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