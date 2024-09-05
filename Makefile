.PHONY: install
install:
	uv sync

.PHONY: clean
clean: clean-build clean-pyc clean-venv

.PHONY: clean-build
clean-build:
	rm -rf build/
	rm -rf testdir/
	rm -rf dist/
	rm -rf .eggs/
	rm -rf src/*.egg-info/
	find . -name '*.egg-info' -delete
	find . -name '*.egg' -delete
	find . -name '*.so' -delete

.PHONY: clean-pyc
clean-pyc:
	find . -name '*.pyc' -delete
	find . -name '*.pyo' -delete
	find . -name '__pycache__' -exec rm -fr {} +

.PHONY: clean-venv
clean-venv:
	rm -f uv.lock
	rm -rf .venv

.PHONY: test
test: install
	uv run --locked python -m unittest discover tests

.PHONY: test_wheel
test_wheel: clean
	@# Build the wheel
	uv build
	@# Now install it in an isolated directory and check that we can import it.
	rm -rf testdir && mkdir testdir
	cd testdir && uv venv
	cd testdir && uv pip install `find ../dist/ -name *.whl`
	testdir/.venv/bin/python -c "from mt2 import mt2"
	rm -r testdir

