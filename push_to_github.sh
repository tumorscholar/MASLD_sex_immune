#!/usr/bin/env bash
# push_to_github.sh -----------------------------------------------------------
# Commit and push this repo to GitHub. Works on the HPC (Apocrita) over HTTPS,
# so no SSH key is required. Repo: tumorscholar/MASLD_sex_immune (branch main).
#
# ONE-TIME SETUP on GitHub:
#   Create a fine-grained Personal Access Token (Settings -> Developer settings
#   -> Personal access tokens) with "Contents: Read and write" on this repo.
#
# USAGE (from the repo root, i.e. the folder that holds this script):
#   export GITHUB_TOKEN=github_pat_xxx            # your token (not stored here)
#   export GITHUB_USER=tumorscholar               # your GitHub username
#   bash push_to_github.sh "your commit message"
#
# The token is used only for this push and is never written to disk or history.
# -----------------------------------------------------------------------------
set -euo pipefail

REPO="tumorscholar/MASLD_sex_immune"
BRANCH="main"
MSG="${1:-Update pipeline and manuscript outputs}"
USER="${GITHUB_USER:-tumorscholar}"

# --- sanity: run from the repo root ------------------------------------------
if [ ! -f "README.md" ] || [ ! -d "scripts" ]; then
  echo "ERROR: run this from the repo root (the folder with README.md and scripts/)." >&2
  exit 1
fi

# --- identity (only sets if unset; harmless to re-run) -----------------------
git config user.name  >/dev/null 2>&1 || git config user.name  "Raju Kumar"
git config user.email >/dev/null 2>&1 || git config user.email "adinraju@gmail.com"

# --- init if this folder is not yet a git repo -------------------------------
if [ ! -d ".git" ]; then
  echo "No .git here - initialising a new repository."
  git init
  git branch -M "$BRANCH"
fi

# --- require the token -------------------------------------------------------
if [ -z "${GITHUB_TOKEN:-}" ]; then
  echo "ERROR: set GITHUB_TOKEN first, e.g.  export GITHUB_TOKEN=github_pat_xxx" >&2
  exit 1
fi

# --- point origin at the authenticated HTTPS URL just for this push ----------
# (token URL is set, used, then reset to the plain URL so the token is not
#  left behind in .git/config)
PLAIN_URL="https://github.com/${REPO}.git"
AUTH_URL="https://${USER}:${GITHUB_TOKEN}@github.com/${REPO}.git"

if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "$AUTH_URL"
else
  git remote add origin "$AUTH_URL"
fi

# --- stage, commit, push -----------------------------------------------------
git add -A
git status --short
git commit -m "$MSG" || echo "(nothing new to commit)"
git branch -M "$BRANCH"

# pull-rebase first if the remote already has commits (safe no-op on a new repo)
git pull --rebase origin "$BRANCH" 2>/dev/null || true

git push -u origin "$BRANCH"

# --- scrub the token back out of .git/config ---------------------------------
git remote set-url origin "$PLAIN_URL"

echo
echo "Done. View at: https://github.com/${REPO}"
echo "(origin reset to the token-free URL: ${PLAIN_URL})"
