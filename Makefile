.PHONY: dist docs

module = pycalculix

dist:
	make dist_docs
	make dist_source

dist_docs:
	sphinx-apidoc -F -A "Justin Black" -o docs pycalculix
	echo 'from pycalculix import __version__' >> docs/conf.py
	echo 'version = __version__' >> docs/conf.py
	echo 'release = version' >> docs/conf.py
	echo "extensions = ['sphinx.ext.autodoc'," >> docs/conf.py
	echo "    'sphinx.ext.napoleon'," >> docs/conf.py
	echo "    'sphinx.ext.todo'," >> docs/conf.py
	echo "    'sphinx.ext.viewcode']" >> docs/conf.py
	sphinx-build -b html docs documentation
	rm -R -f documentation/.doctrees
	zip -r documentation.zip documentation
	mv documentation.zip dist/
	rm -rf docs
	rm -rf documentation

dist_source:
	python3 setup.py sdist --formats=zip
	rm -rf *.egg-info

# add code to clean examples
# add code to make dist_examples
