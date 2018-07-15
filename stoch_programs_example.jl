using StochasticPrograms, Gurobi, LShapedSolvers, Clp


struct SimpleScenario <: AbstractScenarioData
    π::Float64
    d::Vector{Float64}
    q::Vector{Float64}
end


StochasticPrograms.probability(s::SimpleScenario) = s.π

s1 = SimpleScenario(0.4,[500.0,100],[-24.0,-28])
s2 = SimpleScenario(0.6,[300.0,300],[-28.0,-32])


function StochasticPrograms.expected(sds::Vector{SimpleScenario})
    sd = SimpleScenario(1,sum([s.π*s.d for s in sds]),sum([s.π*s.q for s in sds]))
end



sp = StochasticProgram([s1, s2], solver=GurobiSolver())

@first_stage sp = begin
    x = Array{JuMP.Variable}(4)
    x[1] = @variable(model, category=:Int)
    x[2] = @variable(model, category=:Bin)
    x[3] = @variable(model)
    x[4] = @variable(model)

    # @variable(model, x[1:4] >= 3, Int)
    # print(typeof(x))

    # @variable(model, x₂ >= 20)
    @objective(model, Min, sum(x[i] for i in 1:4))
    @constraint(model, x[1] + x[2] >= 120)
    @constraint(model, x[3] + x[4] >= 110)
end

@second_stage sp = begin
    @decision x
    # s = scenario
    # @variable(model, 0 <= y₁ <= s.d[1])
    # @variable(model, 0 <= y₂ <= s.d[2])
    # @objective(model, Min, s.q[1]*y₁ + s.q[2]*y₂)
    # @constraint(model, 6*y₁ + 10*y₂ <= 60*x₁)
    # @constraint(model, 8*y₁ + 5*y₂ <= 80*x₂)
end

# print(sp)


# solve(sp, solver=LShapedSolver(:ls, ClpSolver()))
