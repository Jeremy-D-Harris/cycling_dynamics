#!/bin/bash
# Purpose: Read Comma Separated CSV File
# Author: Vivek Gite under GPL v2.0+
# ------------------------------------------
echo off
sed 1d <Figure3B_Data.csv> Figure3B_Data_fixed.csv
INPUT=Figure3B_Data_fixed.csv
COUNTER=0
OLDIFS=$IFS
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read condition rep1 rep2 rep3

do
  if ((count==0)); then
    #statements

	echo "$condition : $rep1 , $rep2, $rep3"
elif ((count==1)); then
  #statements
  echo "$condition : $rep1 , $rep2, $rep3"
  fi

  COUNTER=$(( COUNTER + 1 ))
	# echo "DOB : $dob"
	# echo "SSN : $ssn"
	# echo "Telephone : $tel"
	# echo "Status : $status"
done < $INPUT
IFS=$OLDIFS
