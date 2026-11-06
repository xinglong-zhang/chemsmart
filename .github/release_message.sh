#!/usr/bin/env bash

set -euo pipefail

CURRENT_TAG="${GITHUB_REF_NAME:-}"

if [[ -z "$CURRENT_TAG" ]]; then
    CURRENT_TAG="$(git describe --tags --exact-match HEAD 2>/dev/null || true)"
fi

if [[ -z "$CURRENT_TAG" ]]; then
    echo "Error: this script must run on a tagged commit." >&2
    exit 1
fi

# Find the closest previous tag in Git history.
PREVIOUS_TAG="$(
    git describe \
        --tags \
        --abbrev=0 \
        "${CURRENT_TAG}^" \
        2>/dev/null || true
)"

echo "## Release ${CURRENT_TAG}"
echo
echo "### Changes"
echo

if [[ -n "$PREVIOUS_TAG" ]]; then
    git log \
        --no-merges \
        --reverse \
        --pretty=format:'* %s' \
        "${PREVIOUS_TAG}..${CURRENT_TAG}"
else
    # First release: include all commits up to this tag.
    git log \
        --no-merges \
        --reverse \
        --pretty=format:'* %s' \
        "${CURRENT_TAG}"
fi

echo
echo
echo "---"
echo "Auto-generated release notes."
