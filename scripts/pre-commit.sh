#!/usr/bin/env bash

# Run unittests as pre-commit hook. Symlink this script
# into your .git/hooks as pre-commit to have tests run
# automatically before committing.
# cd double-beam
# ln -s ../../scripts/pre-commit.sh .git/hooks/pre-commit

echo "Running pre-commit hook"
# save other changes as stash so only the commit is tested
STASH_NAME="pre-commit-$(date +%s)"
git stash save -q --keep-index --include-untracked "$STASH_NAME"
./scripts/run-unittests.sh
UNITTESTS_EXIT_CODE=$?
# unstash
git stash pop -q
if [ $UNITTESTS_EXIT_CODE -ne 0 ]; then
	echo "Tests must pass before commit"
	exit 1
fi
