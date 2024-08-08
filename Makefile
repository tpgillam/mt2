.PHONY: install
install:
	uv sync

.PHONY: clean
clean: clean-build clean-pyc clean-test clean-venv

.PHONY: clean-build
clean-build:
	rm -rf build/
	rm -rf dist/
	rm -rf .eggs/
	find . -name '*.egg-info' -delete
	find . -name '*.egg' -delete

.PHONY: clean-pyc
clean-pyc:
	find . -name '*.pyc' -delete
	find . -name '*.pyo' -delete
	find . -name '__pycache__' -exec rm -fr {} +

.PHONY: clean-test
clean-test:
	rm -rf .pytest_cache

.PHONY: clean-venv
clean-venv:
	rm -f uv.lock
	rm -rf .venv

.PHONY: test
test: install
	uv run --locked pytest

