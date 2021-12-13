#!/usr/bin/env bash
mkdir -p bateria
for n in 2 3 4 10 25 50; do
    for i in $(seq 1 4); do
        mkdir -p bateria/$i
        echo "Calculando a bateria para o exemplo $i, n = $n"
        { time JULIA_LOAD_PATH="./" DRAW=D INSTANCE=$i julia bateria.jl $n | tee bateria/$i/$n.txt; } 2> bateria/$i/$n-time.txt 
    done
done
