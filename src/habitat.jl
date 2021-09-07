struct Habitat{N}
    habs::Vector{Matrix{Float64}}

    function Habitat(habs::Vector{Matrix{Float64}})
        N = length(habs)
        new{N}(habs)
    end
    function Habitat(hab::Matrix{Float64})
        new{1}([hab])
    end
end

length(hab::Habitat{N}) where N = N

getindex(hab::Habitat, i) = getindex(hab.habs, i)

struct HabitatPreference{N, F}
    pref::F

    function HabitatPreference(pref::F) where F<:Tuple
        N = length(pref)
        new{N, F}(pref)
    end
end

function HabitatPreference(prefs...)
    pref = tuple(prefs...)
    HabitatPreference(pref)
end

function HabitatPreference(prefvec::V) where V<:Vector
    pref = tuple(prefvec...)
    HabitatPreference(pref)
end


function (habpref::HabitatPreference{N, F})(hab::Habitat{N}; normalize = true) where {N, F}
    pref = mapreduce(.*, 1:N) do n
        habpref.pref[n].(hab[n])
    end
    if (normalize)
        pref ./= sum(pref)
    end
    pref
end
