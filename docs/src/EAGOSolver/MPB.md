

```@docs
loadproblem!(m::EAGO_NLP_Model, nvar::Int64, ncon::Int64,xL::Vector{Float64},
                xU::Vector{Float64}, gL::Vector{Float64}, gU::Vector{Float64},
                sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)
```

```@docs
loadproblem!(m::EAGO_NLP_Model, nvar::Int64, ncon::Int64,xL::Vector{Float64},
             xU::Vector{Float64}, gL::Vector{Float64}, gU::Vector{Float64},
             f,g)
```

```@docs
    MathProgBase.optimize!(s::EAGO_NLP_Model)
```
