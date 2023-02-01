#!/bin/bash
rm results.csv
# for m in {1..5}
# do
    for n in {1..10}
    do
        ./integrate 0 1 1000000000 $n >> results.csv
    done
# done
