#!/usr/bin/env bash

rm *.trimmed

# combine files from different lanes
cat blue_350bp_L1r1_pe_1.fq blue_350bp_L2_pe_1.fq > blue_350bp_pe_1.fq
cat blue_350bp_L1r2_pe_2.fq blue_350bp_L2_pe_2.fq > blue_350bp_pe_2.fq
cat blue_550bp_L1_pe_1.fq blue_550bp_L2_pe_1.fq > blue_550bp_pe_1.fq
cat blue_550bp_L1_pe_2.fq blue_550bp_L2_pe_2.fq > blue_550bp_pe_2.fq


