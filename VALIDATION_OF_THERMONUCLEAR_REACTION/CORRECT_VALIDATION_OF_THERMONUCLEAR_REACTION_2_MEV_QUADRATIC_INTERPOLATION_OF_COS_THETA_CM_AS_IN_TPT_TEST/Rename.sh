#!/bin/bash

for f in $( find . -type f )
do

  g=$(echo $f | sed 's/T2/T3/')
  mv $f $g

done
