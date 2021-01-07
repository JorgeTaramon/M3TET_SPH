#!/bin/bash
cd /Users/jorge/Dropbox/GEOMAR/m3tet_sph\ \(branch-jorge\)/SETUP_TEST/GPLATES/COASTLINES/16Ma/
for old in *.xy;
do mv $old `basename $old .xy`.txt;
done