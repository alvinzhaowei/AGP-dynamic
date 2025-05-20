# DAGP

**Description**: This is the implementation of DAGP


## Cmake version: 3.28.1


## Parameters for running:
  -graph [address to graph file] -update [address to update file] -wtype [the weights type (listed below)] -a [parameter a] -b [parameter b]

  wtype 0: L-hop transition
  wtype 1: random walk
  wtype 2: geometric distribution
  wtype 3: Possion distribution 

An script example: ./DAGP -graph ./graphs/soc -update ./updates/r_d/soc -wtype 0 -a 0.5 -b 0.5