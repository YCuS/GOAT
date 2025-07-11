#!/usr/bin/env bash
set -e

# Common compiler flags
STD="-std=c++17"
THREAD="-pthread"
# Link in the filesystem library for std::filesystem
FILESYSTEM_LIB="-lstdc++fs"

# Build each tool
g++ $STD $THREAD meme2pwm.cpp           -o meme2pwm
g++ $STD $THREAD get_thr.cpp            -o get_thr
g++ $STD $THREAD search_motif.cpp $FILESYSTEM_LIB -o search_motif

echo "Build successful: meme2pwm, get_thr, search_motif"