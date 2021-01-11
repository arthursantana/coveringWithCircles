#!/usr/bin/env bash
JULIA_LOAD_PATH="./" DRAW=$2 DRAWALL=$3 julia cover.jl $1
