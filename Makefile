.ONESHELL:
ENV_PREFIX=$(shell if conda env list | grep -q chemsmart; then echo "conda run -n chemsmart "; fi)

SHELL := /bin/bash
USE_CONDA ?= true  # Default to true if not explicitly set

.PHONY: help
help:             ## Show the help menu.
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

# === Environment Setup ===

.PHONY: venv
venv:  ## Create a Conda environment if USE_CONDA=true.
	@echo "Debug: USE_CONDA=$(USE_CONDA)"
	@if [ "$(USE_CONDA)" = "true" ]; then \
		echo "Using Conda"; \
		make conda-env; \
	else \
		echo "Using virtualenv"; \
		make virtualenv; \
	fi

.PHONY: conda-env
conda-env:  ## Create a Conda environment.
	@echo "Creating Conda environment 'chemsmart'..."
	@if ! conda env list | grep -q chemsmart; then \
		conda create -n chemsmart python=3.10 -y; \
	fi
	conda run -n chemsmart pip install -U pip
	@echo "Conda environment 'chemsmart' is ready. Activate it with 'conda activate chemsmart'."

.PHONY: virtualenv
virtualenv:  ## Create a virtual environment using virtualenv.
	@if ! command -v python3 > /dev/null; then \
		echo "Python 3 is required but not installed. Exiting."; \
		exit 1; \
	fi
	@if [ ! -d "venv" ]; then \
		python3 -m venv venv; \
	fi
	@source venv/bin/activate && pip install -U pip

# === Project Setup ===

.PHONY: install
install:          ## Install the project in development mode.
	$(ENV_PREFIX)pip install -e .[test]
	$(ENV_PREFIX)pip install types-PyYAML

.PHONY: configure
configure:        ## Run chemsmart configuration interactively.
	@echo "Running chemsmart configuration..."
	chemsmart config
	@echo "Running chemsmart server configuration..."
	chemsmart config server
	@read -p "Enter the path to the Gaussian g16 folder (or press Enter to skip): " gaussian_folder; \
	if [ -n "$$gaussian_folder" ]; then \
		echo "Configuring Gaussian with folder: $$gaussian_folder"; \
		chemsmart config gaussian --folder $$gaussian_folder; \
	else \
		echo "Skipping Gaussian configuration."; \
	fi
	@read -p "Enter the path to the ORCA folder (or press Enter to skip): " orca_folder; \
	if [ -n "$$orca_folder" ]; then \
		echo "Configuring ORCA with folder: $$orca_folder"; \
		chemsmart config orca --folder $$orca_folder; \
	else \
		echo "Skipping ORCA configuration."; \
	fi

.PHONY: show
show: ## Display the current environment information.
	@echo "Current environment:"
	@if [ "$(USE_CONDA)" = "true" ]; then \
		conda env list | grep '*'; \
	fi
	$(ENV_PREFIX)python -V
	$(ENV_PREFIX)python -m site

# === Code Quality ===

.PHONY: fmt
fmt:              ## Format code using black and isort.
	$(ENV_PREFIX)isort --skip pyproject.toml --gitignore .
	$(ENV_PREFIX)black -l 79 chemsmart/
	$(ENV_PREFIX)black -l 79 tests/

.PHONY: lint
lint:             ## Run linters (ruff and black).
	$(ENV_PREFIX)ruff check .
	$(ENV_PREFIX)black -l 79 --check .

# === Testing ===

.PHONY: test
test: lint        ## Run tests and generate coverage report.
	$(ENV_PREFIX)pytest -v --cov-config .coveragerc --cov=chemsmart -l --tb=short --maxfail=1 tests/
	$(ENV_PREFIX)coverage xml
	$(ENV_PREFIX)coverage html

# === Cleanup ===

.PHONY: clean
clean: ## Remove temporary and unnecessary files.
	@find ./ -name '*.pyc' -exec rm -f {} + 2>/dev/null
	@find ./ -name '__pycache__' -exec rm -rf {} + 2>/dev/null
	@find ./ -name 'Thumbs.db' -exec rm -f {} + 2>/dev/null
	@find ./ -name '*~' -exec rm -f {} + 2>/dev/null
	@rm -rf .cache .pytest_cache build dist *.egg-info htmlcov .tox docs/_build 2>/dev/null

