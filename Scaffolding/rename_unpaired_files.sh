#!/usr/bin/env bash

# use python script to rename, then combine

cat unpaired_files | while read line ## first

do

raw=$(echo "$line" | cut -f 1)
renamed=$(echo "$line" | cut -f 2)
pref=$(echo $renamed | cut -f 1 -d '.')

echo $raw
echo $renamed
echo $pref

seq_name_changer.py -i $raw -o $renamed -p $pref\: -k @HWI
 

done


