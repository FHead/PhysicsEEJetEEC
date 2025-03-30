#!/bin/bash

mkdir -p $ProjectBase/Samples/ALEPHMCHungarian

for i in $ProjectBase/Samples/ALEPHMC/*
do
   echo
   echo Processing file \"$i\"...
   echo

	./Execute --Input "${i}" --Output "${i/ALEPHMC/ALEPHMCHungarian}" \
		--Gen tgen --Reco t --Fraction 0.01
done