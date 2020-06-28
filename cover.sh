#!/usr/bin/env bash
JULIA_LOAD_PATH="./" DRAW=$2 julia cover.jl $1
