# Contributing to chemsmart

Welcome to the **chemsmart** project! We appreciate contributions from the community to enhance this Chemistry Simulation and Modeling Automation Toolkit. This guide outlines the steps to contribute effectively.

**Note**: This project requires **Python 3.10**. Instructions are tailored for Linux-based systems (Linux, macOS, BSD, etc.).

## Getting Started

### 1. Fork the Repository
- Navigate to the [chemsmart GitHub repository](https://github.com/xinglong-zhang/chemsmart).
- Click the **Fork** button to create a copy of the repository under your GitHub account.
- Clone your fork to your local machine:
  ```bash
  git clone git@github.com:YOUR_GIT_USERNAME/chemsmart.git
  ```
- Enter the project directory:
  ```bash
  cd chemsmart
  ```
- Add the upstream repository to sync with the main project:
  ```bash
  git remote add upstream https://github.com/xinglong-zhang/chemsmart.git
  ```

### 2. Set Up Your Development Environment
We recommend using a Conda environment to manage dependencies.

- **Create a Conda environment** (named `chemsmart` by default):
  ```bash
  make env
  ```
  If you prefer not to use Conda, you can create a virtual environment:
  ```bash
  make virtualenv
  ```
  **Note**: Conda is strongly recommended for consistency.

- **Activate the environment**:
  ```bash
  conda activate chemsmart
  ```

- **Install the project in development mode**:
  ```bash
  make install-dev
  ```

- **Install pre-commit hooks** (required for developers to enforce code style and quality):
  ```bash
  make pre-commit
  ```

- **Verify the setup** by running tests:
  ```bash
  make test
  ```
  Ensure all tests pass before proceeding.

### 3. Configure chemsmart
- Run the configuration to set up user-specific settings (e.g., paths to Gaussian/ORCA, HPC server details):
  ```bash
  make configure
  ```
- Check the generated `~/.chemsmart` directory to ensure settings match your HPC environment (e.g., SLURM/Torque, scratch directories).
- Source your shell configuration to apply changes:
  ```bash
  source ~/.bashrc
  ```

## Making Contributions

### 4. Create a Branch
- Create a new branch for your contribution:
  ```bash
  git checkout -b my-contribution-branch
  ```
  Use a descriptive branch name (e.g., `feature/add-new-job-type` or `fix/bug-in-parser`).

### 5. Make Your Changes
- Edit files using your preferred editor (PyCharm is recommended).
- If your changes introduce new Python packages, update dependencies:
  ```bash
  make update-deps
  ```

### 6. Format and Lint Your Code
- Code formatting (`black`, `isort`) and linting (`ruff`) are enforced automatically via `pre-commit` hooks when you commit changes, provided you have set up `pre-commit` (see [Set Up Your Development Environment](#2-set-up-your-development-environment)).
- To manually run formatting and linting (optional, if not using `pre-commit`):
  ```bash
  make fmt
  make lint
  ```
- The `pre-commit` hooks also run `make clean` to remove temporary files.

### 7. Test Your Changes
- Write tests for new functionality in the `tests/` directory.
- Run tests to ensure everything works:
  ```bash
  make test
  ```
- Aim for **100% code coverage** in the coverage report. Add tests to your pull request (PR) if coverage is below 100%.

### 8. Document Your Changes
- Update documentation in `docs/` or relevant files if your changes affect usage or functionality.
- Build and preview documentation locally (optional):
  ```bash
  make docs
  ```

### 9. Commit Your Changes
- Use [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) for clear commit messages. Examples:
  ```bash
  git commit -m "feat: add support for new Gaussian job type"
  git commit -m "fix: resolve issue with ORCA input parsing 🎉"
  ```
- Emojis are welcome but optional.
- If `pre-commit` is installed, hooks will automatically format, lint, and clean your code before committing.

### 10. Push and Submit a Pull Request
- Push your branch to your fork:
  ```bash
  git push origin my-contribution-branch
  ```
- On GitHub, create a pull request (PR) from your branch to the `main` branch of `xinglong-zhang/chemsmart`.
- Provide a clear PR description, including the purpose of your changes and any relevant details.
- Wait for CI checks to complete and address any feedback from maintainers.

## Development Tools
The project includes a `Makefile` with useful utilities. View available commands:
```bash
make help
```
Key commands:
- `make env`: Create a Conda environment.
- `make install`: Install the project in development mode.
- `make pre-commit`: Install pre-commit hooks.
- `make fmt`: Format code with `black` and `isort`.
- `make lint`: Run `ruff` linter.
- `make test`: Run tests and generate coverage report.
- `make docs`: Build documentation.
- `make clean`: Remove temporary files.

## Code Style and Quality
- **Formatting**: Use `black` (line length: 79) and `isort` (profile: black).
- **Linting**: Use `ruff` with settings defined in `pyproject.toml`.
- **Testing**: Use `pytest` with 100% coverage goal.
- **Commits**: Follow Conventional Commits for consistency.
- **Pre-commit Hooks**: Enforce style, linting, and cleanup automatically, installed via `make pre-commit`.

## Making a Release (Maintainers Only)
Releases follow [Semantic Versioning](https://semver.org/) (e.g., `X.Y.Z`). To create a release:
1. Ensure all changes are committed and tests pass (`make test`).
2. Run:
   ```bash
   make release
   ```
3. Enter the version number (e.g., `0.1.1`) when prompted.
4. The release process updates the changelog, creates a Git tag, and triggers a GitHub Actions workflow to publish to PyPI.

**Note**: Requires a `PYPI_API_TOKEN` secret in the repository settings, generated from [PyPI](https://pypi.org/account/).

**Caution**: `make release` commits all unstaged changes and modifies the changelog.

## Getting Help
- For issues or feature requests, open an [issue](https://github.com/xinglong-zhang/chemsmart/issues).
- For questions, contact the maintainers via GitHub or email (xinglong.zhang@cuhk.edu.hk).
- Refer to the [README.md](README.md) for setup and usage details.

Thank you for contributing to chemsmart!
