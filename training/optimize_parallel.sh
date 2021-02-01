#!/bin/sh

for i in $(seq 1 10); do
    nohup julia main.jl $i >> errout/$i.log 2>&1 &
done

# To terminate the process,
# $ pgrep -f main.jl | xargs kill -9