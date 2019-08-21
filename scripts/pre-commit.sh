#!/usr/bin/env bash

# Run unittests as pre-commit hook. Symink this script
# into your .git/hooks as pre-commit to have tests run
# automatically before committing.
# cd double-beam
# ln -s ../../scripts/pre-commit.sh .git/hooks/pre-commit

echo "Running pre-commit hook"
./scripts/run-unittests.sh

if [ $? -ne 0 ]; then
	echo "Tests must pass before commit"
	exit 1
fi
