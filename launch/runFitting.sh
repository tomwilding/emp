#!/bin/bash
i=1
for arg in "$@"; do
	i=i+1
	ssh pixel0"$i" << EOF
        cd emp
        R CMD BATCH sources.r
        R CMD BATCH launch/"$arg".r
	&
done