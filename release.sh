#!/bin/bash -e

echo "Building Markdown versions of man pages .."
for SUBCOMMAND in genotype call summarise evolve
do
    echo "Documenting $SUBCOMMAND .."
    cargo run -- $SUBCOMMAND --full-help-roff |pandoc - -t markdown -f man |sed 's/\\\[/[/g; s/\\\]/]/g' |cat <(sed s/SUBCOMMAND/$SUBCOMMAND/ prelude) - >docs/usage/lorikeet-$SUBCOMMAND.md
    echo "Finished documenting $SUBCOMMAND"
done