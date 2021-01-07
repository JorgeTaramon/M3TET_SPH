#!/bin/bash
cd /Users/jorge/Dropbox/GEOMAR/m3tet_sph\ \(branch-jorge\)/SETUP_TEST/GPLATES/COASTLINES/14Ma/
for old in *.xy;
do mv $old `basename $old .xy`.txt;
done