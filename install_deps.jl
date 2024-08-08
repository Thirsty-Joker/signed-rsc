using Pkg
const urls = [
    "https://github.com/danspielman/Laplacians.jl.git#master",
    "https://gitee.com/xhs7700/GeneralGraphs.jl#main",
    "https://gitee.com/xhs7700/GraphDatasets.jl#main",
    "https://gitee.com/xhs7700/LinearAlgebraUtils.jl#main",
]

foreach(url -> Pkg.add(url=url), urls)

Pkg.add([
    "Arpack",
    "CSV",
    "CairoMakie",
    "Combinatorics",
    "DataFrames",
    "DataStructures",
    "GLMakie",
    "Graphs",
    "IJulia",
    "LightGraphs",
    "ProgressMeter",
    "SimpleWeightedGraphs",
    "StatsBase",
    "LinearAlgebra",
    "Random"
])