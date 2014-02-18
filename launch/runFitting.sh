#!/bin/bash
ssh pixel01 << EOF
        cd emp
        R CMD BATCH sources.r
        R CMD BATCH launch/processBlur.r
EOF