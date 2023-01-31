#!/bin/bash
for n in {1..10}
do
    ./integrate 0 1 100000000 $n 
done
