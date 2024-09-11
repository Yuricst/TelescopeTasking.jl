"""Column generation scheme for MTTP sharing passes"""


mutable struct ReducedMasterProblem
    P::Int              # number of schedule patterns
    B::Matrix           # map matrix from targets to patterns, K-by-P
end


mutable struct PricingProblem
    pi::Vector         # dual variables
end