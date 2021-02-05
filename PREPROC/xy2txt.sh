#!/bin/bash
cd /Users/jorge/Tests/SPH_MESH/Trash_00/n37k_nmg3_Vel_BCs_130_100_Myr/
for old in *.xy;
do mv $old `basename $old .xy`.txt;
done