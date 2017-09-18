#!/bin/sh

# The script is for running computations of PID
# There is two types of scripts each set of instances:
#
# GATE instances: XOR, AND, XORAND, RDN, RDNXOR, UNQ, RDNUNQXOR
#
# Gaussian instances: Gaussian distributions with different level of disceretization


# # GATE Instances:

# # Instance Parameters 
# gate="RDNUNQXOR"
# noise=0.001
# sample_size=1000000
# # Solver name 
# solver="ECOS_L"

# # Statistics file prefix
# file_prefix="feas_stats"

# # main
# # Solve the averaged instance 
# julia dot_run_this.jl $solver $gate-$noise-$sample_size.dens $file_prefix"_"$solver"_"$gate"_"$noise"_"$sample_size.csv;

# # Solve the sampled instances 
# for i in $(seq 0 99)
# do
#     julia dot_run_this.jl $solver $gate-$noise-$sample_size.dens-$i $file_prefix"_"$solver"_"$gate"_"$noise"_"$sample_size.csv;
    
# done;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# Gaussian Instances

# Instance Parameters  gauss-1-0.25.dens 
dist="gauss"
lower=1
upper=20
group=5
# Solver name 
solver="ECOS_L"

# Statistics file prefix
file_prefix="feas_stats"

# main

# Solve the sampled instances 
for matrix in $(seq $lower $upper)
do
    for level in  0.25
    do
	julia dot_run_this.jl $solver $dist-$matrix-$level.dens $file_prefix"_"$solver"_"$dist"_"$group.csv;
    done;
done;
