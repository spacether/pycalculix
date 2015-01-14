rmdir /s /q "docs"
sphinx-apidoc -F -e -d 3 -o docs pycalculix
REM include the napoleon extension in conf.py
REM 'sphinxcontrib.napoleon'
cd docs
make html