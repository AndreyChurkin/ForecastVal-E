"""
This function finds the Nucleolus for a coalitional game by solving a sequence of linear programs.
Input: Characteristic function describing the cooperative game.
Output: The Nucleolus solution.

For more details, see: Guajardo, Mario, and Kurt Jörnsten, "Common mistakes in computing the Nucleolus_model_1," European Journal of Operational Research 241, no. 3 (2015).

Andrey Churkin https://andreychurkin.ru/

"""


"""
Note: 
The coalitional structure and the order of coalitions in the characteristic function `V` cannot be arbitrary.
The order must match the structure defined in the `Shapley_matrix` function.
"""


"""
Note: 
Games with positive characteristic functions are cost games.
In such games, players achieve synergy in cost reduction (lower values due to cooperation) and allocate the total cost.

Games with negative characteristic functions are value games.
In such games, players achieve higher values (e.g., profits) when cooperating.

It is important to select the right sign of the characteristic function.
Otherwise, Nucleolus calculations will be incorrect.

Nucleolus.jl uses the cost formulation for cooperative games.
In cost games, the excess is a measure of how satisfied a coalition S is with the cost allocation X. 
The larger the excess of S, the more satisfied the coalition S is.
"""


function Nucleolus(n,V)
    GCvalue = V[n] # The cost of the grand coalition
    u = Shapley_matrix(n) # Get the coalitional structure for this game
    εmain = [] # Optimised excess ε at the iteration
    Xmain = zeros(n) # Allocation chronology

    # Omitting the grand coalition to formulate the coalitional excess constraints:
    V = V[setdiff(1:end, n), :]
    u = u[setdiff(1:end, n), :]

    # # --------------- Nucleolus LP problem #1 (iteration 1) ---------------

    Nucleolus_model_1 = Model()
    @variable(Nucleolus_model_1, ε) # Minimum excess ε among all coalitions
    @variable(Nucleolus_model_1, Xvar[x=1:n]) # Preimputation (cost allocation vector)
    @constraint(Nucleolus_model_1, equ[co=1:2^n-2], ε + sum(u[co,x]*Xvar[x] for x = 1:n) <= V[co])
    @constraint(Nucleolus_model_1, sum(Xvar[x] for x=1:n) == GCvalue)
    @objective(Nucleolus_model_1, Max, ε)

    @suppress begin
        set_optimizer(Nucleolus_model_1, Gurobi.Optimizer)
        optimize!(Nucleolus_model_1)

        println()
        println("Iteration #1")
        println("ε = ", value.(ε))
        println("Xvar = ", value.(Xvar))
        println("Duals = ", shadow_price.(equ))
        println()
    end # @suppress

    FixedConstraintsε = zeros(2^n-2,2^n-2) # The set of all coalitions for which the excess constraint is already satisfied (and is fixed now)
    global vecNonNegativeBuf = trues(2^n-2) # To save fixed values

    for cnstr = 1:2^n-2
        if shadow_price.(equ)[cnstr] >= 0.001
                global vecNonNegativeBuf[cnstr] = false
                global FixedConstraintsε[cnstr] = value.(ε) # These constraints were fixed at interation #1
        end
    end

    global εmain = vcat(εmain,[value.(ε)])
    global Xmain = value.(Xvar')



    # # --------------- The remaining LP problems (iterations) ---------------

    iteration = 2
    fl = 1 # A flag to continue or stop Nucleolus LP problems (iterations)
    while fl >= 1

        global equfalse = []
        global equtrue = []

        for equ=1:2^n-2
            if vecNonNegativeBuf[equ] == true
                global equtrue = vcat(equtrue,[equ])
            else
                global equfalse = vcat(equfalse,[equ])
            end
        end

        Nucleolus_model_2=Model()
        @variable(Nucleolus_model_2, ε2) # minimum excess ε among all coalitions
        @variable(Nucleolus_model_2, Xvar[x=1:n]) # preimputation (cost allocation vector)
        @constraint(Nucleolus_model_2, equ1[co=1:2^n-2; vecNonNegativeBuf[co]==false], FixedConstraintsε[co] + sum(u[co,x]*Xvar[x] for x=1:n) == V[co])
        @constraint(Nucleolus_model_2, equ2[co=1:2^n-2; vecNonNegativeBuf[co]==true], ε2 + sum(u[co,x]*Xvar[x] for x=1:n) <= V[co])

        @constraint(Nucleolus_model_2, sum(Xvar[x] for x=1:n) == GCvalue)
        @objective(Nucleolus_model_2, Max, ε2)

        # println("Solving LP problem #", iteration)

        @suppress begin
            set_optimizer(Nucleolus_model_2,Gurobi.Optimizer)
            optimize!(Nucleolus_model_2)
            # optimize!(Nucleolus_model_2,with_optimizer(Gurobi.Optimizer))

            println()
            println("Iteration #", iteration)
            println("ε = ", value.(ε2))
            println("Xvar = ", value.(Xvar))
            println("equfalse = ", equfalse)
            println("equtrue = ", equtrue)
            println("vecNonNegativeBuf = ", vecNonNegativeBuf)
            println("Duals equ2 = ")
            println(shadow_price.(equ2))
            println()
        end # @suppress
        
        fl = 0
        for cnstr = equtrue
            if shadow_price.(equ2)[cnstr] >= 0.001
                    global vecNonNegativeBuf[cnstr]=false
                    fl = 1
                    global FixedConstraintsε[cnstr]=value.(ε2)
            end
        end

        global εmain = vcat(εmain,[value.(ε2)])
        global Xmain = vcat(Xmain,value.(Xvar'))

        if isapprox(Xmain[iteration,:], Xmain[iteration-1,:]) == true
            fl = 0
        end

        if fl !=0
            iteration += 1
        end
    end

    println()
    println("The Nucleolus calculations are complete. Total number of linear programs solved (iterations) = ", iteration)
    return value.(Xvar')
end






