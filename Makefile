# Thin convenience wrapper for macOS/Linux. Every target delegates to a plain
# command that also works on Windows (see README "Quick start").
# Run inside an activated virtual environment.

PYTHON ?= python

.PHONY: help install format lint test coverage generate-demo-data validate qc run-analysis run-dashboard check

help:
	@echo "Targets: install format lint test coverage generate-demo-data validate qc run-analysis run-dashboard check"

install:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements-dev.txt -e .

format:
	$(PYTHON) -m ruff format src tests
	$(PYTHON) -m ruff check --fix src tests

lint:
	$(PYTHON) -m ruff format --check src tests
	$(PYTHON) -m ruff check src tests

test:
	$(PYTHON) -m pytest

coverage:
	$(PYTHON) -m pytest --cov --cov-report=term-missing

generate-demo-data:
	obw generate-data

validate:
	obw validate

qc:
	obw qc

run-analysis:
	obw run-analysis

run-dashboard:
	obw dashboard

check: lint test
