#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-g++}"
CXXFLAGS=(-std=c++17 -O2 -Wall -Wextra)

"$CXX" "${CXXFLAGS[@]}" meme2pwm.cpp -o meme2pwm
"$CXX" "${CXXFLAGS[@]}" get_thr.cpp -o get_thr

if ! "$CXX" "${CXXFLAGS[@]}" search_motif.cpp -o search_motif; then
  "$CXX" "${CXXFLAGS[@]}" search_motif.cpp -lstdc++fs -o search_motif
fi

echo "Build successful: meme2pwm, get_thr, search_motif"
