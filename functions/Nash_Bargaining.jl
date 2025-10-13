"""
This function finds the Nash bargaining solution for a coalitional game by solving a nonlinear optimisation problem.
Input: Number of players `n`, the coalitional structure of the game and the characteristic function `V`, and the disagreement points `d`.
Output: The Nash bargaining solution.

Andrey Churkin https://andreychurkin.ru/

"""

"""
Note: 
This script can take the full characteristic function (all 2^N coalitions) to derive the feasible set as the Core of the game.
The coalitional structure and the order of coalitions in the characteristic function `V` cannot be arbitrary.
The order must match the structure defined in the `Shapley_matrix` function.
"""

"""
Note: 
This Nash bargaining formulation considers value games and expects positive input values to be allocated.
That is, it is expected that the players achieve some positive profits or savings together.
"""




function Nash_bargaining(n, d; V = nothing, total_value = nothing, objective_formulation = "log")
    println()
    println("Formulating and solving the Nash bargaining problem...")
    if V !== nothing
        total_value = V[n] # Allocating the value of the grand coalition
        println("Using the full coalitional structure and adding the Core constraints.")
    elseif total_value !== nothing
        println("Using the total value only, without the Core constraints.")
    else
        printstyled("❌ WARNING: Provide either `V` or `total_value` to formulate the Nash bargaining problem!"; color = :magenta)
    end

    
    global Nash_bargaining_model = Model()

    @variable(Nash_bargaining_model, Xvar[x=1:n] >= d[x]) #  Payoff allocation, larger than the disagreement points

    # # Efficiency constraint:
    @constraint(Nash_bargaining_model, sum(Xvar[x] for x = 1:n) <= total_value + 10^-6)
    @constraint(Nash_bargaining_model, sum(Xvar[x] for x = 1:n) >= total_value - 10^-6)

    if objective_formulation == "log"
        @NLobjective(Nash_bargaining_model, Max, sum(log(Xvar[i] - d[i]) for i in 1:n))
        println("Log-sum objective function is used for the Nash bargaining model.")
    elseif objective_formulation == "prod"
        @NLobjective(Nash_bargaining_model, Max, prod(Xvar[i] - d[i] for i in 1:n))
        println("Product objective function is used for the Nash bargaining model.")
    else
        printstyled("❌ WARNING: Please select the correct objective function!"; color = :magenta)
    end

    # # Coalitional rationality (the Core constraints):
    if V !== nothing
        u = Shapley_matrix(n) # Get the coalitional structure for this game
        for S = 1:size(u)[1]-1
            if S != n # Do not impose constraints on the grand coalition
                @constraint(Nash_bargaining_model, sum(Xvar[x]*u[S,x] for x = 1:n) >= V[S])
            end
        end
    end

    @suppress begin
        set_optimizer(Nash_bargaining_model, Ipopt.Optimizer)
        set_optimizer_attribute(Nash_bargaining_model, "max_iter", 200) # # To avoid the solver getting stuck for infeasible problems
        optimize!(Nash_bargaining_model)
    end # @suppress

        println("Solver termination status = ", termination_status(Nash_bargaining_model))
        println("Objective value = ",objective_value(Nash_bargaining_model))
        if termination_status(Nash_bargaining_model) == LOCALLY_SOLVED
            printstyled("✅ The Nash bargaining model is solved: the solver converged to an optimal solution."; color = :green)
        else
            printstyled("❌ WARNING: The solver didn't converge! Check the correctness of the problem formulation."; color = :magenta)
        end
        println()
        # println("Xvar = ", value.(Xvar))

    return value.(Xvar'), Nash_bargaining_model
end






