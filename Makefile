.PHONY: help pypi pypi-test clean docs

help:
	@echo "pypi - submit to PyPI server"
	@echo "pypi-test - submit to TestPyPI server"
	@echo "docs - generate Sphinx documentation"
	@echo "clean - remove artifacts"

pypi:
	python setup.py sdist bdist_wheel
	twine upload dist/*

pypi-test:
	python setup.py sdist bdist_wheel
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +
	rm -rf build/
	rm -rf dist/
	rm -rf special.egg-info/
	rm -rf .pytest_cache/
	rm -f .coverage

docs:
	rm -rf docs/api
	python helpers/update_docs_rst_from_README.py
	sphinx-apidoc -o docs special
	cd docs/
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
