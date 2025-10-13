"""
This function checks whether the Core of a cooperative game with `n` players and value function `V` is non-empty.
It solves a feasibility problem (linear programming).

Returns true if the Core is non-empty, false otherwise.

The coalitional structure and the order of coalitions in the characteristic function `V` are defined in the `Shapley_matrix` function.

Andrey Churkin https://andreychurkin.ru/

"""


function Core_nonemptiness(;n, coalitional_structure, V)
    @suppress begin
        global Core_nonemptiness_model = Model(Gurobi.Optimizer)

        @variable(Core_nonemptiness_model, Xvar[x=1:n]) # Payoff allocation to players

        # # Efficiency constraint:
        @constraint(Core_nonemptiness_model, sum(Xvar[x] for x = 1:n) == V[n])

        # # Coalitional rationality (the Core constraints):
        for S = 1:size(coalitional_structure)[1]
            if S != n # Do not impose constraints on the grand coalition
                @constraint(Core_nonemptiness_model, sum(Xvar[x]*coalitional_structure[S,x] for x = 1:n) >= V[S])
            end
        end

        optimize!(Core_nonemptiness_model)
    end # @suppress

        if termination_status(Core_nonemptiness_model) == OPTIMAL
            Core_nonemptiness_result = true
            printstyled("✅ The Core is nonempty."; color = :green)
            println()
        else
            Core_nonemptiness_result = false
            printstyled("❗WARNING: The Core is empty! No feasible allocation exists."; color = :magenta)
            println()
        end    

    return Core_nonemptiness_result
end