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

const _mc_unit_parameters = Dict(
    6 => Dict{Symbol,Vector{ComplexF64}}(
        :alpha => [
            -0.0159969178377436,
            -1.4774047320945315e-10,
            -2.175194571396893,
            -1.435411352267946e-11,
            -0.0047911195018862,
            -1.5738622265225956e-09,
        ],
        :beta => [
            0.7852874768986048im,
            -0.8125036002430518im,
            -1.0787475408258504im,
            -0.6754562463557517im,
            0.8054714129050877im,
        ],
        :w => [
            2.7407079820293235e-05 - 1.114709601324414e-05im,
            -0.4789189669207517 + 0.3988213902999421im,
            6.339100697905341e-06 - 3.531951185814045e-06im,
            0.4816516807336158 - 0.3843678674248487im,
            -1.39688879447589e-06 + 2.446828920026841e-06im,
            0.3828168658192857 - 0.2926144102988278im,
        ],
    ),
    8 => Dict{Symbol,Vector{ComplexF64}}(
        :alpha => [
            -1.0571454550443759e-09,
            -1.6395458268521668e-10,
            -2.696923007439711e-11,
            -2.97808855115383,
            -1.0161736214526636e-09,
            -3.613611453100646e-09,
            -3.533734455676678e-11,
            -3.7308241469738736e-11,
        ],
        :beta => [
            -0.8878551441886222im,
            0.4072216500650863im,
            -0.9963622043839604im,
            -1.4939562087367522im,
            -1.0447733902465357im,
            -0.4549378571694225im,
            0.8475503716040427im,
        ],
        :w => [
            -0.0658317501607036 - 0.2480693025711961im,
            -0.1306190039495692 + 0.0346735860059933im,
            -0.1792101501482174 - 0.6752581452255852im,
            0.0191814455637157 - 0.0050848262355202im,
            0.0977165900714343 + 0.3681497465041551im,
            -0.1357210522537042 + 0.0360130045978403im,
            -0.1063780972151471 - 0.4007719856087376im,
            -0.2912758280346078 + 0.077306047686776im,
        ],
    ),
    10 => Dict{Symbol,Vector{ComplexF64}}(
        :alpha => [
            -0.3430711405230032,
            -8.674886121083773e-05,
            -2.727656090386843,
            -0.7093203697466086,
            -3.2354616763218383e-06,
            -4.5024755928768747e-07,
            -2.785321971202375e-06,
            -9.476253878456016e-05,
            -0.0013707004912771,
            -5.954027014131764e-06,
        ],
        :beta => [
            1.1315098386115716im,
            1.0453507004042732im,
            -1.0753339704488616im,
            0.8349475497916327im,
            -0.6039941768050517im,
            -0.5093377444416272im,
            0.6774846403870297im,
            0.1605800966525998im,
            -0.9507225900537972im,
        ],
        :w => [
            -0.001318528814401 + 0.0004616307985947im,
            0.003320461502283 + 0.0005489766548683im,
            -0.0024046448485234 - 0.0014818661614868im,
            0.0194192645510923 - 0.0355344526619434im,
            -0.0331769757641655 - 0.0120047397905429im,
            0.1037507750668207 - 0.3530710400477936im,
            0.1205615272932349 + 0.0207949826216024im,
            0.1647842706995822 - 0.8174053022264767im,
            -0.1211749374501335 - 0.0044462242353818im,
            0.0471863161980153 - 0.3668399977306861im,
        ],
    ),
)

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
    markovianclosure_parameters(Ω::Real, K::Real, n)

Construct a SemicircleMarkovianClosure object with `n` pseudomodes, asymptotic frequency `Ω`
and coupling `K`, from the universal set of parameters.
"""
function markovianclosure_parameters(Ω::Real, K::Real, n)
    if !haskey(_mc_unit_parameters, n)
        errmsg =
            "Markovian closure not implemented for n = $n. " *
            "Please choose between " *
            join(sort(collect(keys(_mc_unit_parameters))), ", ", " and ") *
            "."
        throw(ArgumentError(errmsg))
    end
    return markovianclosure_parameters(
        Ω,
        K,
        _mc_unit_parameters[n][:alpha],
        _mc_unit_parameters[n][:beta],
        _mc_unit_parameters[n][:w],
    )
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
