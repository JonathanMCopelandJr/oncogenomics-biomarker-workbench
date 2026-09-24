# Thin convenience wrapper for macOS/Linux. Every target delegates to a plain
# command that also works on Windows (see README "Quick start").
# Run inside an activated virtual environment.

PYTHON ?= python

.PHONY: help install format lint test coverage generate-demo-data validate qc run-analysis ml-demo run-dashboard check

help:
	@echo "Targets: install format lint test coverage generate-demo-data validate qc run-analysis ml-demo run-dashboard check"

install:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements-dev.txt -e .

format:
	$(PYTHON) -m ruff format src tests app
	$(PYTHON) -m ruff check --fix src tests app

lint:
	$(PYTHON) -m ruff format --check src tests app
	$(PYTHON) -m ruff check src tests app

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

ml-demo:
	obw ml-demo

run-analysis:
	obw run-analysis

run-dashboard:
	obw dashboard

check: lint test
