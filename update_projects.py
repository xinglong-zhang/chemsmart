import tomlkit

PYPROJECT_PATH = "pyproject.toml"
REQUIREMENTS_PATH = "detected_requirements.txt"

# Load existing pyproject.toml
with open(PYPROJECT_PATH, "r") as f:
    pyproject = tomlkit.parse(f.read())

# Extract current dependencies
current_deps = set(pyproject["project"]["dependencies"])

# Load detected dependencies
with open(REQUIREMENTS_PATH, "r") as f:
    detected_deps = {
        line.strip() for line in f if line.strip() and not line.startswith("#")
    }

# Merge dependencies (avoid duplicates)
updated_deps = sorted(current_deps.union(detected_deps))

# Update pyproject.toml
pyproject["project"]["dependencies"] = updated_deps

# Save the updated file
with open(PYPROJECT_PATH, "w") as f:
    f.write(tomlkit.dumps(pyproject))

print("Updated pyproject.toml with detected dependencies.")
