#!/usr/bin/env bash

open_book=false

# Loop through all arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --open)
            open_book=true
            shift # remove the current argument
            ;;
        *)    # unknown option
            shift # remove the current argument
            ;;
    esac
done

# build the untranslated book
mdbook build
# build the translated german book
MDBOOK_BOOK__LANGUAGE=de mdbook build -d book/de

# Check if the open flag is true
if [ "$open_book" = true ]; then
    # open the book in the browser
    open book/index.html
fi

