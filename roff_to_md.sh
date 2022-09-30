#!/bin/bash -e

echo "Building Markdown versions of man pages .."
for SUBCOMMAND in genotype call summarise consensus
do
    echo "Converting $SUBCOMMAND .."
    sed 's/\\\[/[/g; s/\\\]/]/g' docs/usage/lorikeet-$SUBCOMMAND.wd.md |cat <(sed s/SUBCOMMAND/$SUBCOMMAND/ prelude) - >docs/usage/lorikeet-$SUBCOMMAND.md
    echo "Finished documenting $SUBCOMMAND"

done
rm docs/usage/lorikeet-*.wd.*
sed -i 's/# NAME//' docs/usage/lorikeet-*.md