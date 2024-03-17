# Markovian closure operators

## Constructors

```@docs
SemicircleMarkovianClosure
```

```@docs
markovianclosure_parameters
```

## Properties and parameters

```@docs
length(mc::SemicircleMarkovianClosure)
freqs(mc::SemicircleMarkovianClosure)
innercoups(mc::SemicircleMarkovianClosure)
outercoups(mc::SemicircleMarkovianClosure)
damps(mc::SemicircleMarkovianClosure)
freq(mc::SemicircleMarkovianClosure, j::Int)
innercoup(mc::SemicircleMarkovianClosure, j::Int)
outercoup(mc::SemicircleMarkovianClosure, j::Int)
damp(mc::SemicircleMarkovianClosure, j::Int)
```

## ITensor operators

```@docs
markovianclosure(mc::SemicircleMarkovianClosure, sites::Vector{<:Index}, chain_edge_site::Int)
markovianclosure_adjoint(mc::SemicircleMarkovianClosure, sites::Vector{<:Index}, chain_edge_site::Int, gradefactor::Int)
markovianclosure′
filled_markovianclosure
filled_markovianclosure_adjoint
filled_markovianclosure′
```
