#!/usr/bin/env bash

# use python script to rename, then combine

cat paired_files | while read line

do

raw=$(echo "$line" | cut -f 1)
renamed=$(echo "$line" | cut -f 2)
pref=$(echo $renamed | cut -f 1-4 -d '_')
strand=$(echo "$line" | cut -f 3)

echo $raw
echo $renamed
echo $pref
echo $strand

seq_name_changer.py -i $raw -o $renamed -p $pref -f $strand -k @HWI
 

done


