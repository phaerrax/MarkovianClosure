export SemicircleMarkovianClosure, markovianclosure_parameters
export freq, innercoup, outercoup, damp, freqs, innercoups, outercoups, damps

"""
A SemicircleMarkovianClosure object stores the parameters that define a Markovian closure
which replaces an environment with a semicircle spectral density.
"""
struct SemicircleMarkovianClosure
    frequency
    innercoupling
    outercoupling
    damping
    @doc """
        SemicircleMarkovianClosure(
            ω::Vector{<:Real}, γ::Vector{<:Real}, g::Vector{<:Real}, ζ::Vector{<:Complex}
        )

    SemicircleMarkovianClosure is an aggregate type which stores the parameters needed to
    describe the set of pseudomodes that make up a Markovian closure replacing an
    environment with a semicircle spectral density.
    """
    function SemicircleMarkovianClosure(
        ω::Vector{<:Real}, γ::Vector{<:Real}, g::Vector{<:Real}, ζ::Vector{<:Complex}
    )
        if (length(ω) - 1 != length(g) || length(ω) != length(γ) || length(ω) != length(ζ))
            error("Lengths of input parameters do not match.")
        end
        return new(ω, g, ζ, γ)
    end
end

"""
    markovianclosure_parameters(
        Ω::Real, K::Real,
        α::Matrix{<:Real}, β::Matrix{<:Real}, w::Matrix{<:Real}
    )

Construct a SemicircleMarkovianClosure object with asymptotic frequency `Ω` and coupling
`K`, with the universal parameter set given by `α`, `β` and `w`, each of which is assumed
to be a two-column matrix with the real parts of the parameters in the first column and the
imaginary parts in the second one.
"""
function markovianclosure_parameters(
    Ω::Real, K::Real, α::Matrix{<:Real}, β::Matrix{<:Real}, w::Matrix{<:Real}
)
    return markovianclosure_parameters(
        Ω, K, α[:, 1] .+ im .* α[:, 2], β[:, 1] .+ im .* β[:, 2], w[:, 1] .+ im .* w[:, 2]
    )
end

"""
    markovianclosure_parameters(
        Ω::Real, K::Real,
        α::Vector{<:Complex}, β::Vector{<:Complex}, w::Vector{<:Complex}
    )

Construct a SemicircleMarkovianClosure object with asymptotic frequency `Ω` and coupling
`K`, with the universal parameter set given by `α`, `β` and `w`, which are assumed to be
vectors of complex numbers.
"""
function markovianclosure_parameters(
    Ω::Real, K::Real, α::Vector{<:Complex}, β::Vector{<:Complex}, w::Vector{<:Complex}
)
    frequency = @. Ω - 2K * imag(α)
    damping = @. -4K * real(α)
    innercoupling = @. -2K * imag(β)
    outercoupling = @. K * w
    return SemicircleMarkovianClosure(frequency, damping, innercoupling, outercoupling)
end

"""
    length(mc::SemicircleMarkovianClosure)

Return the length of a Markovian closure, i.e. the number of pseudomodes.
"""
Base.length(mc::SemicircleMarkovianClosure) = Base.length(mc.frequency)

"""
    freqs(mc::SemicircleMarkovianClosure)

Return the list of frequencies of the pseudomodes in the given Markovian closure.
"""
freqs(mc::SemicircleMarkovianClosure) = mc.frequency

"""
    innercoups(mc::SemicircleMarkovianClosure)

Return the list of inner couplings of the pseudomodes in the Markovian closure `mc`,
i.e. the coupling constant between the pseudomodes.
"""
innercoups(mc::SemicircleMarkovianClosure) = mc.innercoupling

"""
    outercoups(mc::SemicircleMarkovianClosure)

Return the list of outer couplings of the pseudomodes in the Markovian closure `mc`,
i.e. the coupling constants with the last point of the TEDOPA chain.
"""
outercoups(mc::SemicircleMarkovianClosure) = mc.outercoupling

"""
    damps(mc::SemicircleMarkovianClosure)

Return the list of damping coefficients of the pseudomodes in the Markovian closure `mc`.
"""
damps(mc::SemicircleMarkovianClosure) = mc.damping

"""
    freq(mc::SemicircleMarkovianClosure, j::Int)

Return the frequency of the `j`-th pseudomode of the Markovian closure `mc`.
"""
freq(mc::SemicircleMarkovianClosure, j::Int) = mc.frequency[j]

"""
    innercoup(mc::SemicircleMarkovianClosure, j::Int)

Return the coupling coefficient between the `j`-th and the `j+1`-th pseudomodes of the
Markovian closure `mc`.
"""
innercoup(mc::SemicircleMarkovianClosure, j::Int) = mc.innercoupling[j]

"""
    outercoup(mc::SemicircleMarkovianClosure, j::Int)

Return the coupling coefficient of the `j`-th pseudomode in the Markovian closure `mc`,
with the last point of the TEDOPA chain.
"""
outercoup(mc::SemicircleMarkovianClosure, j::Int) = mc.outercoupling[j]

"""
    damp(mc::SemicircleMarkovianClosure, j::Int)

Return the damping coefficient of the `j`-th pseudomode in the Markovian closure `mc`.
"""
damp(mc::SemicircleMarkovianClosure, j::Int) = mc.damping[j]
