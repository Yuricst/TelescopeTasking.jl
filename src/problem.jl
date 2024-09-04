"""Struct containing telescope observation scheduling problem information"""


struct TelescopeSchedulingProblem
    n::Int              # number of observation arcs
    m::Int              # number of targets
    A::Matrix           # binary n-by-m matrix associating observation arcs to targets
    T::Matrix           # binary n-by-n matrix denoting feasibility of transition from arc i to j

    function TelescopeSchedulingProblem(n::Int, m::Int, A::Matrix, T::Matrix)
        
        # create model
        model = Model(solver)

        new(n, m, A, T)
    end
end


function Base.show(io::IO, tsp::TelescopeSchedulingProblem)
    println(io, "Telescope Scheduling Problem instance")
    println(io, "    Number of observation arcs: $(tsp.n)")
    println(io, "    Number of targets: $(tsp.m)")
end


"""
Instantiate JuMP model and solve
"""
function solve(problem::TelescopeSchedulingProblem, solver::Gurobi.Optimizer)

    # create model
    model = Model(solver)

    return model
end