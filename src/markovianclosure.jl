export markovianclosure, filled_markovianclosure

using ITensors.SiteTypes: @SiteType_str, SiteType, _sitetypes, commontags

function markovianclosure(
    ::SiteType, ::SemicircleMarkovianClosure, ::Vector{<:Int}, ::Int, ::Int
)
    return nothing
end

"""
    markovianclosure(
        mc::SemicircleMarkovianClosure,
        sites::Vector{<:Index},
        chain_edge_site::Int,
        gradefactor::Int=1,
    )

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`. This closure replaces a chain starting from an empty state.
"""
function markovianclosure(
    mc::SemicircleMarkovianClosure,
    sites::Vector{<:Index},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    @assert length(mc) == length(sites)
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    for st in common_stypes
        ℓ = markovianclosure(st, mc, sitenumber.(sites), chain_edge_site, gradefactor)
        # If the result is something, return that result.
        !isnothing(ℓ) && return ℓ
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"markovianclosure\" function not found for Index " *
            "tags $(tags(sites[1]))",
        ),
    )
end

function markovianclosure(
    st::SiteType"vFermion",
    mc::SemicircleMarkovianClosure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("A†", chain_edge_site, jws..., "A", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("A", chain_edge_site, jws..., "A†", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "A⋅ * ⋅A†"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function markovianclosure(
    st::SiteType"vS=1/2",
    mc::SemicircleMarkovianClosure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, "σ-", site) +
            conj(outercoup(mc, j)) * gkslcommutator("σ-", chain_edge_site, "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        ℓ += gradefactor * damp(mc, j), "σ-⋅ * ⋅σ+", site
        # -0.5 (a† a ρ + ρ a† a)
        ℓ += -0.5damp(mc, j), "N⋅", site
        ℓ += -0.5damp(mc, j), "⋅N", site
    end
    return ℓ
end

function markovianclosure(
    st::SiteType"vElectron",
    mc::SemicircleMarkovianClosure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        # c↑ᵢ† c↑ᵢ₊ₙ = a↑ᵢ† Fᵢ Fᵢ₊₁ ⋯ Fᵢ₊ₙ₋₁ a↑ᵢ₊ₙ
        # c↑ᵢ₊ₙ† c↑ᵢ = -a↑ᵢ Fᵢ Fᵢ₊₁ ⋯ Fᵢ₊ₙ₋₁ a↑ᵢ₊ₙ†
        # c↓ᵢ† c↓ᵢ₊ₙ = a↓ᵢ† Fᵢ₊₁ Fᵢ₊₂ ⋯ Fᵢ₊ₙ a↓ᵢ₊ₙ
        # c↓ᵢ₊ₙ† c↓ᵢ = -a↓ᵢ Fᵢ₊₁ Fᵢ₊₂ ⋯ Fᵢ₊ₙ a↓ᵢ₊ₙ†

        jws = jwstring(; start=chain_edge_site, stop=site)
        # ζⱼ c↑₀† c↑ⱼ (0 = chain edge, j = pseudomode)
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Aup†F", chain_edge_site, jws..., "Aup", site)
        # conj(ζⱼ) c↑ⱼ† c↑₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("AupF", chain_edge_site, jws..., "Aup†", site)
        # ζⱼ c↓₀† c↓ⱼ
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Adn†", chain_edge_site, jws..., "FAdn", site)
        # conj(ζⱼ) c↓ⱼ† c↓₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("Adn", chain_edge_site, jws..., "FAdn†", site)
    end

    # Dissipative part
    for (j, site) in enumerate(sitenumbers)
        # c↑ₖ ρ c↑ₖ† = F₁ ⋯ Fₖ₋₁ a↑ₖ ρ a↑ₖ† Fₖ₋₁ ⋯ F₁
        # Remember that:
        # • Fⱼ = (1 - 2 N↑ₖ) (1 - 2 N↓ₖ);
        # • Fⱼ and aₛ,ₖ commute only on different sites;
        # • {a↓ₖ, a↓ₖ†} = {a↑ₖ, a↑ₖ†} = 1;
        # • Fₖ anticommutes with a↓ₖ, a↓ₖ†, a↑ₖ and a↑ₖ†.
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Aup⋅ * ⋅Aup†"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )
        # c↓ₖ ρ c↓ₖ† = F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ ρ a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "FAdn⋅ * ⋅Adn†F"]
        ℓ += (
            gradefactor * damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...
        )

        # -½ (c↑ₖ† c↑ₖ ρ + ρ c↑ₖ† c↑ₖ) = -½ (a↑ₖ† a↑ₖ ρ + ρ a↑ₖ† a↑ₖ)
        ℓ += -0.5damp(mc, j), "Nup⋅", site
        ℓ += -0.5damp(mc, j), "⋅Nup", site
        # -½ (c↓ₖ† c↓ₖ ρ + ρ c↓ₖ† c↓ₖ) = -½ (a↓ₖ† a↓ₖ ρ + ρ a↓ₖ† a↓ₖ)
        ℓ += -0.5damp(mc, j), "Ndn⋅", site
        ℓ += -0.5damp(mc, j), "⋅Ndn", site
    end
    return ℓ
end

function filled_markovianclosure(
    ::SiteType, ::SemicircleMarkovianClosure, ::Vector{<:Int}, ::Int, ::Int
)
    return nothing
end

"""
    filled_markovianclosure(
        mc::SemicircleMarkovianClosure,
        sites::Vector{<:Index},
        chain_edge_site::Int,
        gradefactor::Int,
    )

Return an OpSum object encoding the Markovian closure operators with parameters given by
`mc`, on sites `sites`, linked to the main TEDOPA/thermofield chain on site
`chain_edge_site`.
This closure replaces a chain starting from a completely filled state.
"""
function filled_markovianclosure(
    mc::SemicircleMarkovianClosure,
    sites::Vector{<:Index},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    @assert length(mc) == length(sites)
    commontags_s = commontags(sites...)
    common_stypes = _sitetypes(commontags_s)
    for st in common_stypes
        ℓ = filled_markovianclosure(
            st, mc, sitenumber.(sites), chain_edge_site, gradefactor
        )
        # If the result is something, return that result.
        !isnothing(ℓ) && return ℓ
        # Otherwise, try again with another type from the initial ones.
    end
    # Return an error if no implementation is found for any type.
    return throw(
        ArgumentError(
            "Overload of \"filled_markovianclosure\" function not found for " *
            "Index tags $(tags(sites[1]))",
        ),
    )
end

function filled_markovianclosure(
    st::SiteType"vFermion",
    mc::SemicircleMarkovianClosure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("A†", chain_edge_site, jws..., "A", site) +
            conj(outercoup(mc, j)) *
            gkslcommutator("A", chain_edge_site, jws..., "A†", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a ρ a†
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "A†⋅ * ⋅A"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # -0.5 (a a† ρ + ρ a a†)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_markovianclosure(
    st::SiteType"vS=1/2",
    mc::SemicircleMarkovianClosure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        ℓ += (
            outercoup(mc, j) * gkslcommutator("σ+", chain_edge_site, "σ-", site) +
            conj(outercoup(mc, j)) * gkslcommutator("σ-", chain_edge_site, "σ+", site)
        )
    end

    for (j, site) in enumerate(sitenumbers)
        # a† ρ a
        ℓ += damp(mc, j), "σ+⋅ * ⋅σ-", site
        # -0.5 (a a† ρ + ρ a a†)
        ℓ += 0.5damp(mc, j), "N⋅", site
        ℓ += 0.5damp(mc, j), "⋅N", site
        ℓ += -damp(mc, j), "Id", site
    end
    return ℓ
end

function filled_markovianclosure(
    st::SiteType"vElectron",
    mc::SemicircleMarkovianClosure,
    sitenumbers::Vector{<:Int},
    chain_edge_site::Int,
    gradefactor::Int=1,
)
    ℓ = spin_chain(st, freqs(mc), innercoups(mc), sitenumbers)

    for (j, site) in enumerate(sitenumbers)
        jws = jwstring(; start=chain_edge_site, stop=site)
        # ζⱼ c↑₀† c↑ⱼ (0 = chain edge, j = pseudomode)
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Aup†F", chain_edge_site, jws..., "Aup", site)
        # conj(ζⱼ) c↑ⱼ† c↑₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("AupF", chain_edge_site, jws..., "Aup†", site)
        # ζⱼ c↓₀† c↓ⱼ
        ℓ +=
            outercoup(mc, j) * gkslcommutator("Adn†", chain_edge_site, jws..., "FAdn", site)
        # conj(ζⱼ) c↓ⱼ† c↓₀
        ℓ +=
            -conj(outercoup(mc, j)) *
            gkslcommutator("Adn", chain_edge_site, jws..., "FAdn†", site)
    end

    # Dissipative part
    for (j, site) in enumerate(sitenumbers)
        # c↑ₖ† ρ c↑ₖ = a↑ₖ† Fₖ₋₁ ⋯ F₁ ρ F₁ ⋯ Fₖ₋₁ a↑ₖ
        # Remember that:
        # • Fⱼ = (1 - 2 N↑ₖ) (1 - 2 N↓ₖ);
        # • Fⱼ and aₛ,ₖ commute only on different sites;
        # • {a↓ₖ, a↓ₖ†} = {a↑ₖ, a↑ₖ†} = 1;
        # • Fₖ anticommutes with a↓ₖ, a↓ₖ†, a↑ₖ and a↑ₖ†.
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Aup†⋅ * ⋅Aup"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)
        # c↓ₖ† ρ c↓ₖ = a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁ ρ F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ
        opstring = [repeat(["F⋅ * ⋅F"], site - 1); "Adn†F⋅ * ⋅FAdn"]
        ℓ += (damp(mc, j), collect(Iterators.flatten(zip(opstring, 1:site)))...)

        # c↑ₖ c↑ₖ† = F₁ ⋯ Fₖ₋₁ a↑ₖ a↑ₖ† Fₖ₋₁ ⋯ F₁ =
        #          = F₁² ⋯ Fₖ₋₁² a↑ₖ a↑ₖ† =
        #          = a↑ₖ a↑ₖ† =
        #          = 1 - N↑ₖ
        ℓ += -damp(mc, j), "Id", site
        ℓ += 0.5damp(mc, j), "Nup⋅", site
        ℓ += 0.5damp(mc, j), "⋅Nup", site
        # c↓ₖ c↓ₖ† = F₁ ⋯ Fₖ₋₁ Fₖa↓ₖ a↓ₖ†Fₖ Fₖ₋₁ ⋯ F₁ =
        #          = F₁² ⋯ Fₖ₋₁² Fₖa↓ₖ a↓ₖ†Fₖ =
        #          = Fₖa↓ₖ a↓ₖ†Fₖ =
        #          = Fₖ² a↓ₖ a↓ₖ† =
        #          = 1 - N↓ₖ
        ℓ += -damp(mc, j), "Id", site
        ℓ += 0.5damp(mc, j), "Ndn⋅", site
        ℓ += 0.5damp(mc, j), "⋅Ndn", site
    end
    return ℓ
end
