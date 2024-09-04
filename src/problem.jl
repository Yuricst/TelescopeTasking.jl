"""Struct containing telescope observation scheduling problem information"""


struct TelescopeSchedulingProblem
    n::Int                      # number of observation arcs
    m::Int                      # number of targets
    A::Matrix                   # binary n-by-m matrix associating observation arcs to targets
    T::Matrix                   # binary n-by-n matrix for transition feasibility from arc i to j
    X::Array{VariableRef, 1}    # binary m-vector denoting whether target k is observed
    Y::Array{VariableRef, 1}    # binary n-vector denoting whether observation arc i is selected
    model::Model                # JuMP model

    function TelescopeSchedulingProblem(
        passes::Vector{VisiblePass}, num_exposure::Int, slew_rate::Number, solver)
        # construct target allocation matrix
        designators = unique([pass.tle.international_designator for pass in passes])
        m = length(designators)
        n = length(passes)
        A = zeros(Int, n, m)
        for (i, pass) in enumerate(passes)
            j = findfirst(x -> x == pass.tle.international_designator, designators)
            A[i,j] = 1
        end

        # construct transition feasibility matrix
        T = zeros(Int, n, n)
        for i in 1:n
            for j in 1:n
                if passes[i].tf < passes[j].t0      # FIXME - need to incorporate slew rate
                    T[i,j] = 1
                end
            end
        end

        # create model
        model = Model(solver)
        @variable(model, Y[1:n], Bin);      # whether observation arc i is selected
        @variable(model, X[1:m], Bin);      # whether target k is observed (sufficiently many times)
        
        @constraint(model, sufficient_exposure[k=1:m], 
                    sum(A[i,k] * Y[i] for i in 1:n) >= num_exposure * X[k])

        @constraint(model, transition_feasibility[i=1:n-1, j=i+1:n], 
                    Y[i] + Y[j] <= 1 + T[i,j])

        @objctive(model, Max, sum(X))
        new(n, m, A, T, X, Y, model)
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
function solve(problem::TelescopeSchedulingProblem)


    return model
end