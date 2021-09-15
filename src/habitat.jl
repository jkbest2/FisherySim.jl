struct Habitat{T, N}
    # Store habitats as a Tuple because there probably won't be very many
    # habitats in most cases, and it allows disparate element types (e.g. some
    # Float64 and some Bool).
    habs::T

    function Habitat(habs::T) where {T<:Tuple}
        N = length(habs)
        chksz = size(habs[1])
        for i in 2:N
            _check_habsize(habs[i], chksz) || error("Habitat $i size does not match")
        end
        new{T, N}(habs)
    end
end

function _check_habsize(hab, chksz)
    sz = size(hab)
    all(sz .== chksz)
end

function Habitat(habs::Vector{<:Matrix})
    Habitat(tuple(habs...))
end

function Habitat(habs...)
    Habitat(tuple(habs...))
end

length(hab::Habitat{T, N}) where {T, N} = N
size(hab::Habitat) = size(hab.habs[1])

getindex(hab::Habitat, i) = getindex(hab.habs, i)

struct HabitatPreference{F, N}
    pref::F

    function HabitatPreference(pref::F) where F<:Tuple
        N = length(pref)
        new{F, N}(pref)
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

function (habpref::HabitatPreference{F, N} where F)(hab::Habitat{T, N} where T; normalize = true) where {N}
    pref = mapreduce(.*, 1:N) do n
        habpref.pref[n].(hab[n])
    end
    if (normalize)
        pref ./= sum(pref)
    end
    pref
end
