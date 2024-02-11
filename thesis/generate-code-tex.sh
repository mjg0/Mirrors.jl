#!/bin/bash

# Generate a tex file ($2) containing the contents of each file listed in $1, with Julia syntax highlighting

cat "$1" | while read line; do
    echo "\\section{$line}\\label{$(basename "$line")}"
    echo '\begin{minted}{julia}'
    cat "../$line"
    echo '\end{minted}'
done > "$2"