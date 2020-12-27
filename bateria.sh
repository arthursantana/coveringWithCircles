#!/usr/bin/env bash
mkdir -p bateria
for n in $(seq 2 25); do
    JULIA_LOAD_PATH="./" DRAW=D julia bateria.jl $n > bateria/$n.txt
done
