#!/usr/bin/python

import sys
sys.path.append("py_mods")
import filter_hmmscan
filter_hmmscan.tab2gff("COGs/orf_tab/Pgab.tab")
