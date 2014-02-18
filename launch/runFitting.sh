#!/bin/bash
i=0
for arg in "$@"; do
        let "i+=1"
        echo ssh pixel0"$i" 'cd emp;R CMD BATCH sources.r;R CMD BATCH launch/'"$arg"'.r'
        ssh pixel0"$i" 'cd emp;R CMD BATCH sources.r;R CMD BATCH launch/'"$arg"'.r' &
done