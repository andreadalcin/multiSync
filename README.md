# Contents

## Testing details
Tested on MATLAB R2020b

## Dependencies
The code for MPLS is taken from https://github.com/yunpeng-shi/MPLS and adapted to match our notation

## Benchmarks
demo_multisyncSO3.m
Runs a benchmark to evaluate the performance of Multisync w.r.t. edge averaging for rotation synchronization on SO(3).


demo_partitioned.m
Runs a benchmark to evaluate the performance of partitioned synchronization using both Multisync and edge averaging on 1DSfM data sets. 
Data sets can be downloaded from https://www.cs.cornell.edu/projects/1dsfm/ and should be added to the 'datasets_snavely' folder. Set ops.data_name to the name of the data set at line 24.
