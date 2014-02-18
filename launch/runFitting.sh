#!/bin/bash
i=0
if [ $# = 0 ]; then
	echo Not Enough Args
fi
for arg in "$@"; do
        let "i+=1"
        ssh pixel0"$i" 'cd emp;R CMD BATCH -q launch/'"$arg"'.r output/'"$arg"'.txt' &
done