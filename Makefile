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
	find . -path ./.venv -prune -o -name '*.egg-info' -exec rm -rf {} +
	find . -path ./.venv -prune -o -name '*.egg' -exec rm -rf {} +
	find . -path ./.venv -prune -o -name '*.so' -exec rm -rf {} +

.PHONY: clean-pyc
clean-pyc:
	find . -path ./.venv -prune -o -name '*.pyc' -exec rm -rf {} +
	find . -path ./.venv -prune -o -name '*.pyo' -exec rm -rf {} +
	find . -path ./.venv -prune -o -name '__pycache__' -exec rm -rf {} +

.PHONY: clean-venv
clean-venv:
	rm -rf .venv

.PHONY: test
test: install
	uv run --locked python -m unittest discover tests

.PHONY: lint
lint: install
	-uv run --locked ruff check --fix
	-uv run --locked ruff format

.PHONY: typecheck
typecheck: install
	uv run --locked pyright

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

