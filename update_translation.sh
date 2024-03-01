#!/usr/bin/env bash

# extract source strings
MDBOOK_OUTPUT='{"xgettext": {"pot-file": "messages.pot", "granularity": 1, "depth": 1}}' mdbook build -d po
# update translation files
msgmerge --update po/de.po po/messages.pot

