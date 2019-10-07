module FisherySim

using Distributions
using Rtweedie
using StatsBase
using StatsFuns         # For `logit` and `logistic`; worth it?
using LinearAlgebra
using Random
using PDMats

import Base: rand, +, -, step, sum, getindex, setindex!, size, length, copy
import Rtweedie: Tweedie
import StatsBase: sample, cov

include("fisherydomain.jl")
export AbstractFisheryDomain,
       DiscreteFisheryDomain,
       GriddedFisheryDomain,
       size,
       length,
       sample
       ## Don't export:
       ## calculate_distances,
       ## map_symm

include("covkernels.jl")
export AbstractCovarianceKernel,
       ExpCov,
       Matérn32Cov,
       Matern32Cov,
       AR1,
       cov

include("bathymetry.jl")
export BathymetryModel,
       Bathymetry,
       rand

include("pop_dynamics.jl")
export PopulationDynamicsModel,
       PopState,
       Schaefer,
       PellaTomlinson,
       StochasticProduction,
       vecstate,
       step,
       sum,
       setindex!,
       copy

include("movement.jl")
export MovementModel, eqdist, approx_eqdist

include("vessels.jl")
export AbstractTargetingBehavior,
       RandomTargeting,
       PreferentialTargeting,
       target,
       Catchability,
       Vessel,
       Catch,
       logistic,
       CPUE,
       +, -,
       fish!,
       getindex,
       Fleet,
       vessels

include("simulation.jl")
export simulate

end # module
