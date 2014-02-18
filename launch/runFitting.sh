#!/bin/bash
i=0
for arg in "$@"; do
        let "i+=1"
        ssh pixel0"$i" 'cd emp;R CMD BATCH launch/'"$arg"'.r' &
done