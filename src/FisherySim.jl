module FisherySim

using Distributions
using TweedieDistributions
using StatsBase
using LinearAlgebra
using Random
using PDMats
using Arpack
using NeutralLandscapes

import Base:
    rand,
    +, -, *,
    step,
    sum,
    getindex, setindex!,
    size, length, eachindex,
    copy
import StatsBase:
    sample,
    cov
import Distributions:
    location

include("fisherydomain.jl")
export
    AbstractFisheryDomain,
    DiscreteFisheryDomain,
    GriddedFisheryDomain,
    size,
    length,
    eachindex,
    sample

include("covkernels.jl")
export
    AbstractCovarianceKernel,
    ExpCov,
    Matérn32Cov,
    Matern32Cov,
    AR1,
    cov

# include("matrixlognormal.jl")
# export
#     MatrixLogNormal,
#     location

include("domaindistributions.jl")
export
    AbstractDomainDistribution,
    DomainDistribution,
    MultiDomainDistribution,
    BlendedDomainDistribution,
    ClassifiedDomainDistribution,
    domain,
    getindex,
    length

include("habitat.jl")
export
    Habitat,
    getindex,
    length,
    HabitatPreference

include("bathymetry.jl")
export
    BathymetryModel,
    Bathymetry,
    rand

include("pop_dynamics.jl")
export
    PopulationDynamicsModel,
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
export
    MovementModel,
    eqdist,
    MovementRate

include("targeting.jl")
export
    AbstractTargetingBehavior,
    RandomTargeting,
    FixedTargeting,
    AbstractPreferentialTargeting,
    PreferentialTargeting,
    DynamicPreferentialTargeting,
    target,
    reset!

include("catchability.jl")
export
    AbstractCatchability,
    Catchability,
    DensityDependentCatchability,
    *

include("vessels.jl")
export
    Vessel,
    Catch,
    CPUE,
    +, -,
    fish!,
    getindex,
    Fleet,
    vessels

include("simulation.jl")
export
    simulate

end # module
