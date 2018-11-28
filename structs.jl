struct Subproblem
    model::JuMP.Model
    constr_π::Vector{JuMP.ConstraintRef}
    constr_λ::Vector{JuMP.ConstraintRef}
    y::JuMP.JuMPArray{JuMP.Variable,1,Tuple{Vector{String}}}
end

struct StageData
    cols::Vector{Vector{String}}
    rhs::Dict{String,Float64}
    obj::Dict{SubString{String},Float64}
    rows::Vector{Any}
    constrs::Dict{SubString{String},Dict{String,Float64}}
    L::Float64
end

struct Stoch_Scenario
    # p probability, d right side, q obj, a coeffs, comps either "E"(quality), "G"(reater or equal) or "L"(ess or equal)
    p::Float64
    h::Vector{Float64}
    q::Vector{Float64}
    T:: Matrix{Float64}
    W:: Matrix{Float64}
    comps::Vector{String}
end
