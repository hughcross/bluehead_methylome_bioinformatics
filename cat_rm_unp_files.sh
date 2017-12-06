#!/usr/bin/env bash

rm *.unpaired

# combine files from different lanes
cat blue_350bp_L1_unp_1.fq blue_350bp_L2_unp_1.fq blue_350bp_L1_unp_2.fq blue_350bp_L2_unp_2.fq > blue_350bp_unpaired.fq
cat blue_550bp_L1_unp_1.fq blue_550bp_L2_unp_1.fq blue_550bp_L1_unp_2.fq blue_550bp_L2_unp_2.fq > blue_550bp_unpaired.fq
