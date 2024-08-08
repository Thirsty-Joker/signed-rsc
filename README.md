# Signed Subgraph Centrality

### Description
We implement three algorithms for resolvent subgraph centrality, including SolverRSC, ForestRSC-S, ForestRSC-E

### Directory Structure
- algorithm.jl: SolverRSC, ForestRSC-S, ForestRSC-E 
- main.jl: entry point
- graph.jl: read graphs
- new_graph.jl: construct a new graph
- install_deps.jl: install all dependencies

### Datasets
Datasets are from Koblenz Network Collection and SNAP.

### Run in command
```bash
cd signed-rsc
julia install_deps.jl
julia main.jl
```
