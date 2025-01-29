# How to develop on this project

chemsmart welcomes contributions from the community.

**You need PYTHON 3.10!**

This instructions are for linux base systems. (Linux, MacOS, BSD, etc.)
## Setting up your own fork of this repo.

- On github interface click on `Fork` button.
- Clone your fork of this repo. `git clone git@github.com:YOUR_GIT_USERNAME/chemsmart.git`
- Enter the directory `cd chemsmart`
- Add upstream repo `git remote add upstream https://github.com/xinglong-zhang/chemsmart`

## Setting up your own conda environment

Follow the Installation Guide on README. Briefly, run
```bash
make env
```
to create a conda environment.

## Install the project in develop mode

Then run
```bash
make install
```
to install the project in develop mode.

## Run the tests to ensure everything is working

Run 
```bash
make test
```
to run the tests.

## Create a new branch to work on your contribution

Run
```bash
git checkout -b my_contribution_branch
```

## Make your changes

Edit the files using your preferred editor (we recommend PyCharm IDE).

## Format the code

Run
Run
```bash
make fmt
```
to format the code. We use `isort` to sort the import statements alphabetically and `black` formatter.

## Run the linter

Run
```bash
make lint
```
to run the linter. We use `ruff` as linter to check for potential bugs, syntax errors, and stylistic issues.

## Test your changes

Please include relevant tests for the new functionalities that you are introducing. Then, run
```bash
make test
```
to run these tests and make sure that they pass.

Ensure code coverage report shows `100%` coverage, add tests to your PR.

## Build the docs locally (optional)

Run 
```bash
make docs
```
to build the docs.

Ensure your new changes are documented.

## Commit your changes

This project uses [conventional git commit messages](https://www.conventionalcommits.org/en/v1.0.0/).

Example: `fix(package): update setup.py arguments 🎉` (emojis are fine too)

## Push your changes to your fork

Run 
```bash
git push origin my_contribution_branch
```

## Submit a pull request

On github interface, click on `Pull Request` button.

Wait CI to run and one of the developers will review your PR.

## Makefile utilities

This project comes with a `Makefile` that contains a number of useful utility.

```bash 
❯ make
Usage: make <target>

Targets:
help            Show the help menu.
env             Create a Conda environment if USE_CONDA=true.
conda-env       Create a Conda environment.
virtualenv      Create a virtual environment using virtualenv.
install         Install the project in development mode.
configure       Run chemsmart configuration interactively.
show            Display the current environment information.
fmt             Format code using black and isort.
lint            Run linters (ruff).
test            Run tests and generate coverage report.
clean           Remove temporary and unnecessary files.
```

## Making a new release (TO BE IMPLEMENTED)

This project uses [semantic versioning](https://semver.org/) and tags releases with `X.Y.Z`
Every time a new tag is created and pushed to the remote repo, github actions will
automatically create a new release on github and trigger a release on PyPI.

For this to work you need to setup a secret called `PIPY_API_TOKEN` on the project settings>secrets, 
this token can be generated on [pypi.org](https://pypi.org/account/).

To trigger a new release all you need to do is:

1. If you have changes to add to the repo
    * Make your changes following the steps described above.
    * Commit your changes following the [conventional git commit messages](https://www.conventionalcommits.org/en/v1.0.0/).
2. Run the tests to ensure everything is working.
4. Run `make release` to create a new tag and push it to the remote repo.

the `make release` will ask you the version number to create the tag, ex: type `0.1.1` when you are asked.

> **CAUTION**:  The make release will change local changelog files and commit all the unstaged changes you have.
