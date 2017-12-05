module FisherySim

using Distributions
using StatsBase
using Distances

include("bathymetry.jl")
export BathymetryModel, Bathymetry, rand
end # module
