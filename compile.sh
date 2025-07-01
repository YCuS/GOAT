#!/bin/bash
set -e

# Common compiler flags
STD="-std=c++17"
THREAD="-pthread"

# 编译时加上 -pthread
g++ $STD $THREAD meme2pwm.cpp   -o meme2pwm
g++ $STD $THREAD get_thr.cpp    -o get_thr
g++ $STD $THREAD search_motif.cpp -o search_motif

echo "Build successful: meme2pwm, get_thr, search_motif"