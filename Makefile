# Detect the operating system
OS := $(shell uname -s 2>/dev/null || echo Windows)
ifeq ($(OS),Windows)
    SHELL := cmd
    ENV_PREFIX := $(if $(shell where conda >nul 2>&1 && conda env list | findstr chemsmart >nul 2>&1),conda activate chemsmart && ,)
    SEP := \\
    RM := del /Q
    RMDIR := rmdir /S /Q
    ECHO := echo
    NULL := nul
else
    SHELL := /bin/bash
    ENV_PREFIX := $(shell if conda env list | grep -q chemsmart; then echo "conda run -n chemsmart "; fi)
    SEP := /
    RM := rm -f
    RMDIR := rm -rf
    ECHO := echo
    NULL := /dev/null
endif

USE_CONDA ?= true  # Default to true if not explicitly set
MAKEFILE_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
CHEMSMART_PATH := $(MAKEFILE_DIR)chemsmart$(SEP)cli$(SEP)chemsmart  # Use platform-specific separator

# === Help messages for make ===

.PHONY: help
ifeq ($(OS),Windows)
help:             ## Show the help menu.
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@type $(MAKEFILE_LIST) | findstr /R "^[a-zA-Z_-]*:.*## " | for /F "tokens=1,2 delims=##" %%a in ('more') do @echo %%a                    %%b
else
help:             ## Show the help menu.
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
endif

# === Environment Setup ===

.PHONY: env
env:  ## Create a Conda environment if USE_CONDA=true.
	@echo Debug: USE_CONDA=$(USE_CONDA)
ifeq ($(OS),Windows)
	@if "$(USE_CONDA)"=="true" ( \
		$(ECHO) "Using Conda" && $(MAKE) conda-env \
	) else ( \
		$(ECHO) "Using virtualenv" && $(MAKE) virtualenv \
	)
else
	@if [ "$(USE_CONDA)" = "true" ]; then \
		$(ECHO) "Using Conda"; \
		$(MAKE) conda-env; \
	else \
		$(ECHO) "Using virtualenv"; \
		$(MAKE) virtualenv; \
	fi
endif

.PHONY: conda-env
conda-env:  ## Create or update the Conda environment using environment.yml.
	@echo Managing Conda environment 'chemsmart' with environment.yml...
ifeq ($(OS),Windows)
	@if not exist environment.yml ( \
		$(ECHO) "Error: environment.yml not found in $(MAKEFILE_DIR). Please create it first." && exit 1 \
	)
	@conda env list | findstr chemsmart >$(NULL) && ( \
		$(ECHO) "Updating existing 'chemsmart' environment..." && \
		conda env update -n chemsmart -f environment.yml --prune \
	) || ( \
		$(ECHO) "Creating new 'chemsmart' environment..." && \
		conda env create -f environment.yml \
	)
else
	@if [ ! -f environment.yml ]; then \
		$(ECHO) "Error: environment.yml not found in $(MAKEFILE_DIR). Please create it first."; \
		exit 1; \
	fi
	@if conda env list | grep -q chemsmart; then \
		$(ECHO) "Updating existing 'chemsmart' environment..."; \
		conda env update -n chemsmart -f environment.yml --prune; \
	else \
		$(ECHO) "Creating new 'chemsmart' environment..."; \
		conda env create -f environment.yml; \
	fi
endif
	@echo Conda environment 'chemsmart' is ready. Activate it with 'conda activate chemsmart'.

.PHONY: virtualenv
virtualenv:  ## Create a virtual environment using virtualenv.
ifeq ($(OS),Windows)
	@where python3 >$(NULL) 2>&1 || ( $(ECHO) "Python 3 is required but not installed. Exiting." && exit 1 )
	@if not exist "venv" ( python3 -m venv venv )
	@call venv\Scripts\activate.bat && pip install -U pip
else
	@if ! command -v python3 >$(NULL); then \
		$(ECHO) "Python 3 is required but not installed. Exiting."; \
		exit 1; \
	fi
	@if [ ! -d "venv" ]; then \
		python3 -m venv venv; \
	fi
	@. venv/bin/activate && pip install -U pip
endif

# === Project Setup ===

.PHONY: install
install:          ## Install the project in development mode.
	$(ENV_PREFIX)pip install .  # Normal users (runtime only)
	$(ENV_PREFIX)pip install types-PyYAML

.PHONY: install-dev
install:          ## Install the project in development mode.
	$(ENV_PREFIX)pip install -e .[dev,test]  # install dependencies in dev and test too
	$(ENV_PREFIX)pip install types-PyYAML

.PHONY: pre-commit
pre-commit:       ## Install pre-commit hooks to enforce code style and quality.
	$(ENV_PREFIX)pre-commit install
	@echo Pre-commit hooks installed. They will run automatically on each commit.

.PHONY: configure
configure:        ## Run chemsmart configuration interactively.
ifeq ($(OS),Windows)
	@echo Running chemsmart configuration...
	$(ENV_PREFIX)python $(CHEMSMART_PATH) config
	@echo Running chemsmart server configuration...
	$(ENV_PREFIX)python $(CHEMSMART_PATH) config server || ( $(ECHO) "Error: chemsmart server configuration failed." && exit 1 )
	@echo "Interactive path prompts are not implemented for Windows shell in this Makefile."
	@echo "Please run these manually if needed:"
	@echo "  chemsmart config gaussian --folder <path>"
	@echo "  chemsmart config orca --folder <path>"
	@echo "  chemsmart config nciplot --folder <path>"
else
	@echo Running chemsmart configuration...
	$(ENV_PREFIX)python $(CHEMSMART_PATH) config
	@echo Running chemsmart server configuration...
	$(ENV_PREFIX)python $(CHEMSMART_PATH) config server || ( $(ECHO) "Error: chemsmart server configuration failed." && exit 1 )
	@read -p "Enter the path to the Gaussian g16 folder (or press Enter to skip): " gaussian_folder; \
	if [ -n "$$gaussian_folder" ]; then \
		$(ECHO) "Configuring Gaussian with folder: $$gaussian_folder"; \
		$(ENV_PREFIX)python $(CHEMSMART_PATH) config gaussian --folder "$$gaussian_folder"; \
	else \
		$(ECHO) "Skipping Gaussian configuration."; \
	fi; \
	read -p "Enter the path to the ORCA folder (or press Enter to skip): " orca_folder; \
	if [ -n "$$orca_folder" ]; then \
		$(ECHO) "Configuring ORCA with folder: $$orca_folder"; \
		$(ENV_PREFIX)python $(CHEMSMART_PATH) config orca --folder "$$orca_folder"; \
	else \
		$(ECHO) "Skipping ORCA configuration."; \
	fi; \
	read -p "Enter the path to the NCIPLOT folder (or press Enter to skip): " nciplot_folder; \
	if [ -n "$$nciplot_folder" ]; then \
		$(ECHO) "Configuring NCIPLOT with folder: $$nciplot_folder"; \
		$(ENV_PREFIX)python $(CHEMSMART_PATH) config nciplot --folder "$$nciplot_folder"; \
	else \
		$(ECHO) "Skipping NCIPLOT configuration."; \
	fi
endif

.PHONY: show
show: ## Display the current environment information.
	@echo Current environment:
ifeq ($(OS),Windows)
	@if "$(USE_CONDA)"=="true" conda env list ^| findstr "*"
else
	@if [ "$(USE_CONDA)" = "true" ]; then conda env list | grep '*'; fi
endif
	$(ENV_PREFIX)python -V
	$(ENV_PREFIX)python -m site

# === Code Quality ===

.PHONY: update-deps
update-deps:          ## Automatically update new packages that are added in the codes
	@echo Updating additional dependencies to pyproject.toml file...
	$(ENV_PREFIX)python $(CHEMSMART_PATH) update deps
	@echo Reinstalling chemsmart package...
	$(ENV_PREFIX)pip install -e .[test]

.PHONY: fmt
fmt:              ## Format code using black and isort.
	$(ENV_PREFIX)isort --skip pyproject.toml --gitignore .
	$(ENV_PREFIX)black -l 79 .

.PHONY: lint
lint:             ## Run linters (ruff).
	$(ENV_PREFIX)ruff check . --fix

# === Testing ===

.PHONY: coverage-clean
coverage-clean:   ## Remove any stale coverage files prior to running tests.
ifeq ($(OS),Windows)
	-@for /R . %%f in (.coverage*) do @$(RM) "%%f" 2>$(NULL)
else
	-@rm -f .coverage .coverage.* 2>/dev/null
endif

.PHONY: test
test: lint coverage-clean ## Run tests and generate coverage report (robust to corrupt shards).
	$(ENV_PREFIX)pytest -v --cov-config .coveragerc --cov=chemsmart --cov-branch -l --tb=short --maxfail=1 tests/
	-$(ENV_PREFIX)coverage combine  # Portable error ignoring: - so a bad shard cannot fail the job 
	-$(ENV_PREFIX)coverage xml
	-$(ENV_PREFIX)coverage html

# === Docs ===
.PHONY: docs-clean docs-lint docs-fmt docs

docs-clean: ## Clean documentation artifacts.
	$(MAKE) -C docs clean

docs-lint: ## Lint docs with doc8 and Sphinx (warnings as errors)
	@echo "==> Installing doc8, sphinx, and all doc extensions..."
	$(ENV_PREFIX)pip install -q doc8 sphinx myst_parser sphinx_rtd_theme sphinx-autodoc-typehints
	@echo "==> Running doc8..."
	$(ENV_PREFIX)doc8 --max-line-length=120 --ignore-path docs/build docs/source
	@echo "==> Running sphinx-build -W..."
	$(ENV_PREFIX)sphinx-build -W -b html docs/source docs/build

# Format all .rst files in docs/source using rstfmt
docs-fmt: ## Auto-format reStructuredText with rstfmt.
	@echo "==> Installing rstfmt (if needed)..."
	$(ENV_PREFIX)pip install -q rstfmt
	@echo "==> Running rstfmt..."
	# Format recursively; --in-place edits files
	$(ENV_PREFIX)rstfmt -w 120 docs/source

docs: ## Build documentation (HTML).
	$(MAKE) -C docs html


# === Cleanup ===

.PHONY: clean
clean: ## Remove temporary and unnecessary files.
ifeq ($(OS),Windows)
	@for /R . %%f in (*.pyc) do @$(RM) "%%f" 2>$(NULL)
	@for /D /R . %%d in (__pycache__) do @if exist "%%d" $(RMDIR) "%%d" 2>$(NULL)
	@for /R . %%f in (Thumbs.db) do @$(RM) "%%f" 2>$(NULL)
	@for /R . %%f in (*~) do @$(RM) "%%f" 2>$(NULL)
	@$(RMDIR) .cache .pytest_cache build dist *.egg-info htmlcov .tox .coverage.* docs\_build 2>$(NULL)
else
	@find ./ -name '*.pyc' -exec rm -f {} + 2>/dev/null
	@find ./ -name '__pycache__' -exec rm -rf {} + 2>/dev/null
	@find ./ -name 'Thumbs.db' -exec rm -f {} + 2>/dev/null
	@find ./ -name '*~' -exec rm -f {} + 2>/dev/null
	@rm -rf .cache .pytest_cache build dist *.egg-info htmlcov .tox .coverage.* docs/_build 2>/dev/null
endif

