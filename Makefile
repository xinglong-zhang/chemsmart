.ONESHELL:
ENV_PREFIX=$(shell conda run -n chemsmart which python 2>/dev/null || echo "not_found")
CONDA_ENV_NAME=chemsmart

.PHONY: help
help:             ## Show the help.
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@fgrep "##" Makefile | fgrep -v fgrep

.PHONY: show
show:             ## Show the current environment.
	@echo "Current environment:"
	@if [ "$(ENV_PREFIX)" != "not_found" ]; then \
		echo "Using Conda environment '$(CONDA_ENV_NAME)'"; \
		$(ENV_PREFIX) --version; \
		$(ENV_PREFIX) -m site; \
	else \
		echo "Conda environment '$(CONDA_ENV_NAME)' is not active or doesn't exist."; \
	fi

.PHONY: install
install:          ## Install the project in dev mode.
	@echo "Installing the project in editable mode..."
	@$(ENV_PREFIX) -m pip install -e .[test]

.PHONY: fmt
fmt:              ## Format code using black & isort.
	$(ENV_PREFIX) -m isort chemsmart/
	$(ENV_PREFIX) -m black -l 79 chemsmart/
	$(ENV_PREFIX) -m black -l 79 tests/

.PHONY: lint
lint:             ## Run pep8, black, mypy linters.
	$(ENV_PREFIX) -m flake8 chemsmart/
	$(ENV_PREFIX) -m black -l 79 --check chemsmart/
	$(ENV_PREFIX) -m black -l 79 --check tests/
	$(ENV_PREFIX) -m mypy --ignore-missing-imports chemsmart/

.PHONY: test
test: lint        ## Run tests and generate coverage report.
	$(ENV_PREFIX) -m pytest -v --cov-config .coveragerc --cov=chemsmart -l --tb=short --maxfail=1 tests/
	$(ENV_PREFIX) -m coverage xml
	$(ENV_PREFIX) -m coverage html

.PHONY: watch
watch:            ## Run tests on every change.
	ls **/**.py | entr $(ENV_PREFIX) -m pytest -s -vvv -l --tb=long --maxfail=1 tests/

.PHONY: clean
clean:            ## Clean unused files.
	@find ./ -name '*.pyc' -exec rm -f {} \;
	@find ./ -name '__pycache__' -exec rm -rf {} \;
	@find ./ -name 'Thumbs.db' -exec rm -f {} \;
	@find ./ -name '*~' -exec rm -f {} \;
	@rm -rf .cache
	@rm -rf .pytest_cache
	@rm -rf .mypy_cache
	@rm -rf build
	@rm -rf dist
	@rm -rf *.egg-info
	@rm -rf htmlcov
	@rm -rf .tox/
	@rm -rf docs/_build

.PHONY: virtualenv
virtualenv:       ## Create or update the Conda environment from environment.yml.
	@echo "Creating or updating the Conda environment '$(CONDA_ENV_NAME)' from environment.yml..."
	conda env create --name $(CONDA_ENV_NAME) --file environment.yml || \
	conda env update --name $(CONDA_ENV_NAME) --file environment.yml --prune
	@echo
	@echo "Environment '$(CONDA_ENV_NAME)' is ready."
	@echo "Activate it with: conda activate $(CONDA_ENV_NAME)"
	@echo "Run 'conda activate $(CONDA_ENV_NAME)' to use the environment."

.PHONY: release
release:          ## Create a new tag for release.
	@echo "WARNING: This operation will create a version tag and push to GitHub"
	@read -p "Version? (provide the next x.y.z semver) : " TAG
	@echo "$${TAG}" > chemsmart/VERSION
	@$(ENV_PREFIX) -m gitchangelog > HISTORY.md
	@git add chemsmart/VERSION HISTORY.md
	@git commit -m "release: version $${TAG} 🚀"
	@echo "Creating git tag : $${TAG}"
	@git tag $${TAG}
	@git push -u origin HEAD --tags
	@echo "GitHub Actions will detect the new tag and release the new version."

.PHONY: docs
docs:             ## Build the documentation.
	@echo "Building documentation..."
	@$(ENV_PREFIX) -m mkdocs build
	URL="site/index.html"; xdg-open $$URL || sensible-browser $$URL || x-www-browser $$URL || gnome-open $$URL || open $$URL

.PHONY: init
init:             ## Initialize the project based on an application template.
	@./.github/init.sh

