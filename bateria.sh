#!/usr/bin/env bash
mkdir -p bateria
#for n in 2 3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 75 100; do
#for n in 2 3 4 5 6 7 8 9 10 15 20 25; do
#for n in 30 35 40 45; do
#for n in 50 75 100; do
for n in 75; do
    #for i in $(seq 1 4); do
    for i in 3; do
        mkdir -p bateria/$i
        echo "Calculando a bateria para o exemplo $i, n = $n"
        { time JULIA_LOAD_PATH="./" DRAW=D INSTANCE=$i julia bateria.jl $n | tee bateria/$i/$n.txt; } 2> bateria/$i/$n-time.txt 
    done
done
