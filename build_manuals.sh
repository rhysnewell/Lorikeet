#!/bin/bash -e

echo "Building ROFF versions of man pages .."
for SUBCOMMAND in genotype call summarise consensus
do
    echo "Documenting $SUBCOMMAND .."
    cargo run -- $SUBCOMMAND --full-help-roff > docs/usage/lorikeet-$SUBCOMMAND.wd.roff
    echo "Finished documenting $SUBCOMMAND"
done