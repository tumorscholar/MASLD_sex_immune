#!/usr/bin/env bash
# push_to_github.sh -----------------------------------------------------------
# Commit and push this repo to GitHub over SSH.
# Remote: git@github.com:tumorscholar/MASLD_sex_immune.git  (branch main)
#
# Run from the repo root, on any machine whose SSH key is on your GitHub account
# (your laptop or the HPC). Usage:
#     bash push_to_github.sh "your commit message"
# -----------------------------------------------------------------------------
set -e
MSG="${1:-Update pipeline: end-to-end runner, tables, in-house Fig 7 provenance}"

# make sure we are at the repo root
if [ ! -d ".git" ]; then
  echo "Run this from the repo root (the folder that contains .git)." >&2
  exit 1
fi

# make sure the remote is the SSH one (safe to re-run)
git remote set-url origin git@github.com:tumorscholar/MASLD_sex_immune.git 2>/dev/null || \
  git remote add origin git@github.com:tumorscholar/MASLD_sex_immune.git

git add -A
git status --short
git commit -m "$MSG" || echo "nothing to commit"
git branch -M main
git push -u origin main

echo
echo "Done. View at: https://github.com/tumorscholar/MASLD_sex_immune"
