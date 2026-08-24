#!/bin/bash

julia -O3 -C native --math-mode=fast -g0 --project=../../.. -t4 sim.jl $1