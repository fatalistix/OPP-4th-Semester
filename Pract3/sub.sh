#!/bin/sh

echo "=============================================="
echo "parallel with 1 proc:"
echo
./parallel_new 1
./parallel_new 1
./parallel_new 1
echo "=============================================="
echo "parallel with 2 proc:"
echo
./parallel_new 2
./parallel_new 2
./parallel_new 2
echo "=============================================="
echo "parallel with 4 proc:"
echo
./parallel_new 4
./parallel_new 4
./parallel_new 4
echo "=============================================="
echo "parallel with 8 proc:"
echo
./parallel_new 8
./parallel_new 8
./parallel_new 8
echo "=============================================="
echo "parallel with 16 proc:"
echo
./parallel_new 16
./parallel_new 16
./parallel_new 16
echo "=============================================="
echo "parallel with 24 proc:"
echo
./parallel_new 24
./parallel_new 24
./parallel_new 24
echo "=============================================="
