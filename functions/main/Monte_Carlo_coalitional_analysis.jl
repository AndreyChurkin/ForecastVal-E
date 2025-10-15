"""
This script performs a Monte Carlo analysis to assess how the accuracy of forecasts affects operational outcomes and the forecast valuations.
For each iteration, the total avoided cost (savings) due to forecast sharing is calculated, and the forecast valuation is performed.

Several solution concepts are implemented to allocate the savings among forecast providers, including:
- The Shapley value
- Nucleolus
- The Core
- Leave-one-out approach

Andrey Churkin https://andreychurkin.ru/

"""


cd(dirname(@__FILE__))
pwd()

using PowerModels, JuMP, Ipopt, Gurobi
using Suppressor
using JLD2


# # Select the solver to use in the OPF formulations:
optimize_with = "Gurobi"
# optimize_with = "Ipopt"

# # Select the .m network data file to solve OPF:

# file = "../../cases/case33bw (pure).m" # Note: the original case33bw has parsing problems
# file = "../../cases/case33bw (mod) v1.0.m" # <-- voltage at node 1 fixed to 1.0 pu
file = "../../cases/case33bw (mod) v1.1.m" # <-- case33bw used in the paper "Tracing, Ranking and Valuation of Aggregated DER Flexibility in Active Distribution Networks"




network_data = parse_file(file)

L = length(network_data["branch"])
N = length(network_data["bus"])

total_pd = sum(item["pd"] for item in values(network_data["load"]))
total_qd = sum(item["qd"] for item in values(network_data["load"]))

println()
println("Total P load of the network: ",round(total_pd*network_data["baseMVA"],digits=3)," MW")
println("Total Q load of the network: ",round(total_qd*network_data["baseMVA"],digits=3)," MVAr")

# # check the radial structure of the network:
for i = 1:N
      global in_sum = 0 # number of flows "to" the bus
      for l = 1:length(network_data["branch"])
            if network_data["branch"][string(l)]["t_bus"] == i
                  if network_data["branch"][string(l)]["br_status"] == 1
                        in_sum += 1
                  end
            end
      end
      if in_sum >= 2
            println()
            printstyled("Warning: the network might be not radial!"; color = :magenta)
            println()
      end
end

# # check the number of generators (sources)：
if length(network_data["gen"]) == 1
      global gen_bus = network_data["gen"]["1"]["gen_bus"] # this is the bus with the main generator (TSO/DSO interface)
else
      println()
      printstyled("Warning: there is more than one generator (source) in the network!"; color = :magenta)
      println()
end

# # check the lines connected to the source:
global gen_line = 0
global count_gen_line = 0
for l = 1:length(network_data["branch"])
      if network_data["branch"][string(l)]["f_bus"] == gen_bus
            if network_data["branch"][string(l)]["br_status"] == 1
                  global gen_line = l # this is the line connected to the main generator (TSO/DSO interface)
                  global count_gen_line += 1
            end
      end
      if network_data["branch"][string(l)]["t_bus"] == gen_bus
            if network_data["branch"][string(l)]["br_status"] == 1
                  global gen_line = l # this is the line connected to the main generator (TSO/DSO interface)
                  global count_gen_line += 1
            end
      end
end
if count_gen_line >= 2
      println()
      printstyled("Warning: there is more than one line connected to the generator (source)!"; color = :magenta)
      println()
end



# # Select whether to impose MVA ratings on branches:
""" Note: some .m cases have "rate_a"=0 for lines, thus, imposing line limits will cause infeasibility (or parsing errors)"""
# include_line_limits = true
include_line_limits = false



# # Select whether to fix the voltage (at 1.0 pu) at the primary substation:
primary_substation_voltage_fixed = true
# primary_substation_voltage_fixed = false



# # Data of power supply from the source bus, e.g., from the transmission system (TSO):
TSO_data = Dict("bus_i" => gen_bus,
                "cost_DA" => 100, # <-- [EUR/MWh]
                "cost_BM_up" => 300, # <-- [EUR/MWh]
                "cost_BM_down" => 200 # <-- [EUR/MWh]
)



# # Use this part to introduce flexible generators:

# Creating a dictionary for flexible units' data:
flex_data = Dict("bus" =>
            Dict(string(i) =>
                    Dict("bus_i" => i, "flex_p_up" => 0.0, "flex_p_down" => 0.0, "flex_q_up" => 0.0, "flex_q_down" => 0.0, 
                    "cost_DA_up" => 0.0, "cost_DA_down" => 0.0, "cost_BM_up" => 0.0, "cost_BM_down" => 0.0)
                    for i = 1:N ) 
)

# # # Defining flexible units in the network:
flex_data["bus"]["7"]["flex_p_up"] = 1/network_data["baseMVA"]
flex_data["bus"]["7"]["flex_p_down"] = 1/network_data["baseMVA"]
flex_data["bus"]["7"]["flex_q_up"] = 0/network_data["baseMVA"] # Note: we focus only on active power regulation and costs
flex_data["bus"]["7"]["flex_q_down"] = 0/network_data["baseMVA"]
flex_data["bus"]["7"]["cost_DA_up"] = 200 # <-- [EUR/MWh]
flex_data["bus"]["7"]["cost_DA_down"] = 200 # <-- [EUR/MWh]
flex_data["bus"]["7"]["cost_BM_up"] = 400 # <-- [EUR/MWh]
flex_data["bus"]["7"]["cost_BM_down"] = 400 # <-- [EUR/MWh]




# # Use this part to define EVA and their forecasts:

# creating a dictionary for EVA data:
EVA_data = Dict("bus" =>
            Dict(string(i) =>
                    Dict("bus_i" => i, 
                        "load_forecast" => 0.0, 
                        "load_observation" => 0.0,
                        "status_quo_forecast" => 0.0
                        )
                    for i = 1:N ) 
)

"""
Note: status quo forecast is the best forecast initially available to the DSO.
This is the load forecast used if no additional forecast is provided by EVA.
We need such a forecast to perform data valuation (to estimate the difference between with and without EVA forecast).
The term status quo is not perfect and may be changed in the future.
"""

# # Defining EVAs in the network:

EVA_data["bus"]["22"]["load_forecast"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["22"]["load_observation"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["22"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

EVA_data["bus"]["32"]["load_forecast"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["32"]["load_observation"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["32"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

EVA_data["bus"]["14"]["load_forecast"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["14"]["load_observation"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["14"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

EVA_data_copy = deepcopy(EVA_data) # <-- we will need this for iterative load modifications

# # Define the players of a cooperative game: names of EVAs (their buses), and their order:

# Players = [14] # <-- one player
# Players = [32 14] # <-- 2 players
Players = [22 32 14] # <-- 3 players



# # How many Monte Carlo simulations to perform:
Monte_Carlo_runs = 10 # (was 1000 in the paper)

"""
Note: we will only vary load_observation and load_forecast in each simulation
The status_quo_forecast and the value of all loads are kept intact
"""


# # Generate new random simulations or read existing saved Monte Carlo simulations:
# read_predefined_Monte_Carlo_data = true
read_predefined_Monte_Carlo_data = false


""" 
If true, the load forecasts and observations will be read from a file. If false, new random values will be generated. 
"""

if read_predefined_Monte_Carlo_data == true
    global predefined_Monte_Carlo_data = load("../../results/Monte_Carlo_coalitional_analysis_results_1000_positive_runs.jld2")
    # global Monte_Carlo_runs = length(predefined_Monte_Carlo_data["Monte_Carlo_results_for_export"])
end



# # ------------ Data valuation starts here ------------

# # Shapley values for each iteration will be saved here:
Sh_iterative_results = zeros(Monte_Carlo_runs, length(Players))

EVA_forecasts_iterative_results = zeros(0,3)
EVA_observations_iterative_results = zeros(0,3)


# # Functions for defining the coalitional structure and calculating the Shapley value:
include("../valuation/Shapley_matrix.jl")
include("../valuation/Shapley_value.jl")
include("../valuation/allocation_rescaling.jl")
# # Function to check whether the Core of a cooperative game is non-empty (returns true if non-empty, false otherwise):
include("../valuation/Core_nonemptiness_check.jl")
# # Function to calculate the Nucleolus:
include("../valuation/Nucleolus.jl")

u = Shapley_matrix(length(Players)) # get the coalitional structure

u_buses = zeros(size(u,1), N) # relation between buses and coalitions: EVAs at which buses belong to which coalitions

for coalition = 1:2^length(Players)
      for player = 1:length(Players)
            if u[coalition,player] == 1
                  global u_buses[coalition,Players[player]] = 1
            end
      end
end

# # Writing down dispatch results for each iteration (EV load value):
iterative_results_for_export = Dict{Int, Dict}()

@time begin

for Monte_Carlo_iteration = 1:Monte_Carlo_runs
    skip_this_iteration = false

    println()
    println("Monte Carlo simulation #",Monte_Carlo_iteration)

    V_iteration = zeros(2^length(Players)) # The characteristic function for this iteration

    # # Using random load forecasts and observations:
    EVA_forecasts = []
    EVA_observations = []

    if read_predefined_Monte_Carlo_data == false # No Monte Carlo data given, generate new
        for player_i in Players
            # # # ±50% load forecasts:
            # EVA_data["bus"][string(player_i)]["load_forecast"] = (rand() + 0.5)*EVA_data_copy["bus"][string(player_i)]["load_forecast"]
            # EVA_forecasts = vcat(EVA_forecasts, EVA_data["bus"][string(player_i)]["load_forecast"]*network_data["baseMVA"])

            # # # ±25% load forecasts:
            EVA_data["bus"][string(player_i)]["load_forecast"] = (rand()/2 + 0.75)*EVA_data_copy["bus"][string(player_i)]["load_forecast"]
            EVA_forecasts = vcat(EVA_forecasts, EVA_data["bus"][string(player_i)]["load_forecast"]*network_data["baseMVA"])


            # # # ±50% load observations:
            # EVA_data["bus"][string(player_i)]["load_observation"] = (rand() + 0.5)*EVA_data_copy["bus"][string(player_i)]["load_observation"]
            # EVA_observations = vcat(EVA_observations, EVA_data["bus"][string(player_i)]["load_observation"]*network_data["baseMVA"])
        
            # # # Fixed load observations:
            EVA_data["bus"][string(player_i)]["load_observation"] = 1.0*EVA_data_copy["bus"][string(player_i)]["load_observation"]
            EVA_observations = vcat(EVA_observations, EVA_data["bus"][string(player_i)]["load_observation"]*network_data["baseMVA"])
        end
    else # Read Monte Carlo data
        if !haskey(predefined_Monte_Carlo_data["Monte_Carlo_results_for_export"][Monte_Carlo_iteration], :coalitional_analysis)
            for (i,player_i) in enumerate(Players)
                EVA_data["bus"][string(player_i)]["load_forecast"] = predefined_Monte_Carlo_data["Monte_Carlo_results_for_export"][Monte_Carlo_iteration][:EVA_forecasts][i]/network_data["baseMVA"]
                EVA_forecasts = vcat(EVA_forecasts, EVA_data["bus"][string(player_i)]["load_forecast"]*network_data["baseMVA"])

                EVA_data["bus"][string(player_i)]["load_observation"] = predefined_Monte_Carlo_data["Monte_Carlo_results_for_export"][Monte_Carlo_iteration][:EVA_observations][i]/network_data["baseMVA"]
                EVA_observations = vcat(EVA_observations, EVA_data["bus"][string(player_i)]["load_observation"]*network_data["baseMVA"])
            end
        else
            skip_this_iteration = true # Skip - coalitional analysis already performed for this iteration.
        end
    end


    if skip_this_iteration == false
        global EVA_forecasts_iterative_results = vcat(EVA_forecasts_iterative_results, EVA_forecasts')
        global EVA_observations_iterative_results = vcat(EVA_observations_iterative_results, EVA_observations')

        # # First, we need to estimate the cost with no EVAs forecasts (only with DSO's status quo forecast):

        # println()
        # println("Solving the day-ahead and balancing dispatches for the case with no EVA forecasts")

        @suppress begin
            include("../optimisation/DA_dispatch_without_forecasts.jl") # <-- day-ahead dispatch formulation
            optimize!(OPFmodel_DA)
        end # @suppress

        global DA_dispatch_solution = Dict(
            "gen" => Dict(),
            "flex" => Dict(),
            "total_DA_cost" => Dict()
            )

        global DA_dispatch_solution["gen"] = Dict("p_opt_DA" => value.(pg[1]), "q_opt_DA" => value.(qg[1]))

        for bus = 1:N
        DA_dispatch_solution["flex"][string(bus)] = 
            Dict("pg_flex_up_opt_DA" => value.(pg_flex_up)[bus],
                    "pg_flex_down_opt_DA" => value.(pg_flex_down)[bus],                                   
        )
        end

        global DA_dispatch_solution["total_DA_cost"] = Dict("EUR/h" => objective_value(OPFmodel_DA)*network_data["baseMVA"])

        @suppress begin
            include("../optimisation/BM_dispatch.jl") # <-- balancing dispatch formulation
            optimize!(OPFmodel_BM)
        end # @suppress

        global BM_dispatch_solution = Dict(
        "gen" => Dict(),
        "flex" => Dict(),
        "total_DA_cost" => Dict()
        )

        global BM_dispatch_solution["gen"] = Dict("p_opt_BM_up" => value.(pg_BM_up[1]), "q_opt_BM_up" => value.(pg_BM_up[1]),
                                    "p_opt_BM_down" => value.(pg_BM_down[1]), "q_opt_BM_down" => value.(pg_BM_down[1])        
        )

        for bus = 1:N
        BM_dispatch_solution["flex"][string(bus)] = 
            Dict("pg_flex_up_opt_BM" => value.(pg_flex_BM_up)[bus],
                    "pg_flex_down_opt_BM" => value.(pg_flex_BM_down)[bus]                                   
        )
        end

        global BM_dispatch_solution["total_BM_cost"] = Dict("EUR/h" => objective_value(OPFmodel_BM)*network_data["baseMVA"])

        # # We will need this to compute the avoided cost due to EVA forecasts:
        global total_DA_and_BM_cost_without_forecast = DA_dispatch_solution["total_DA_cost"]["EUR/h"] + BM_dispatch_solution["total_BM_cost"]["EUR/h"]

        # println("Total cost (DA + BM) = ",total_DA_and_BM_cost_without_forecast, " EUR/h")
        # println()


        for cl = 1:2^length(Players)
            global coalition = cl

            # println("Coalition # ",coalition)

            # # First, the day-ahead dispatch is formulated and solved:

            @suppress begin
                    include("../optimisation/DA_dispatch_for_coalition.jl") # <-- day-ahead dispatch formulation
                    optimize!(OPFmodel_DA)
            end # @suppress
                
            global DA_dispatch_solution = Dict(
            "gen" => Dict(),
            "flex" => Dict(),
            "total_DA_cost" => Dict()
            )

            global DA_dispatch_solution["gen"] = Dict("p_opt_DA" => value.(pg[1]), "q_opt_DA" => value.(qg[1]))

            for bus = 1:N
            DA_dispatch_solution["flex"][string(bus)] = 
                    Dict("pg_flex_up_opt_DA" => value.(pg_flex_up)[bus],
                        "pg_flex_down_opt_DA" => value.(pg_flex_down)[bus],                                   
            )
            end

            global DA_dispatch_solution["total_DA_cost"] = Dict("EUR/h" => objective_value(OPFmodel_DA)*network_data["baseMVA"])

            if cl == length(Players)
                global GC_DA_dispatch_cost = objective_value(OPFmodel_DA)*network_data["baseMVA"]

                global GC_DA_voltages = sqrt.(value.(w))

                global GC_pg_DA = value.(pg)*network_data["baseMVA"]
                global GC_qg_DA = value.(qg)*network_data["baseMVA"]

                global GC_pg_flex_DA_up = value.(pg_flex_up)*network_data["baseMVA"]
                global GC_pg_flex_DA_down = value.(pg_flex_down)*network_data["baseMVA"]
                global GC_qg_flex_DA_up = value.(qg_flex_up)*network_data["baseMVA"]
                global GC_qg_flex_DA_down = value.(qg_flex_down)*network_data["baseMVA"]

            end


            # # Then, the balancing market dispatch is formulated and solved:

            @suppress begin
                    include("../optimisation/BM_dispatch.jl") # <-- balancing dispatch formulation
                    optimize!(OPFmodel_BM)
            end # @suppress

            global BM_dispatch_solution = Dict(
            "gen" => Dict(),
            "flex" => Dict(),
            "total_DA_cost" => Dict()
            )

            global BM_dispatch_solution["gen"] = Dict("p_opt_BM_up" => value.(pg_BM_up[1]), "q_opt_BM_up" => value.(pg_BM_up[1]),
                                            "p_opt_BM_down" => value.(pg_BM_down[1]), "q_opt_BM_down" => value.(pg_BM_down[1])        
            )

            for bus = 1:N
            BM_dispatch_solution["flex"][string(bus)] = 
                    Dict("pg_flex_up_opt_BM" => value.(pg_flex_BM_up)[bus],
                        "pg_flex_down_opt_BM" => value.(pg_flex_BM_down)[bus]                                   
            )
            end

            global BM_dispatch_solution["total_BM_cost"] = Dict("EUR/h" => objective_value(OPFmodel_BM)*network_data["baseMVA"])

            if cl == length(Players)
                global GC_BM_dispatch_cost = objective_value(OPFmodel_BM)*network_data["baseMVA"]

                global GC_BM_voltages = sqrt.(value.(w))

                global GC_pg_flex_BM_up = value.(pg_flex_BM_up)*network_data["baseMVA"]
                global GC_pg_flex_BM_down = value.(pg_flex_BM_down)*network_data["baseMVA"]
                global GC_qg_flex_BM_up = value.(qg_flex_BM_up)*network_data["baseMVA"]
                global GC_qg_flex_BM_down = value.(qg_flex_BM_down)*network_data["baseMVA"]

                global GC_pg_BM_up = value.(pg_BM_up)*network_data["baseMVA"]
                global GC_pg_BM_down = value.(pg_BM_down)*network_data["baseMVA"]
                global GC_qg_BM_up = value.(qg_BM_up)*network_data["baseMVA"]
                global GC_qg_BM_down = value.(qg_BM_down)*network_data["baseMVA"]
            end

            # println("Total cost (DA + BM) = ",DA_dispatch_solution["total_DA_cost"]["EUR/h"] + BM_dispatch_solution["total_BM_cost"]["EUR/h"], " EUR/h")
            
            global coalition_avoided_cost = total_DA_and_BM_cost_without_forecast - (DA_dispatch_solution["total_DA_cost"]["EUR/h"] + BM_dispatch_solution["total_BM_cost"]["EUR/h"])
            
            # println("Avoided cost due to the coalition's aggregate forecast = ", coalition_avoided_cost)
            # println()
            
            global V_iteration[coalition] = coalition_avoided_cost

        end

        GC_avoided_cost = V_iteration[length(Players)] # avoided cost of the grand coalition

        # println()
        println("Avoided cost of the Grand Coalition = ", GC_avoided_cost)

        Core_nonemptiness_check = Core_nonemptiness(n = length(Players), coalitional_structure = u, V = V_iteration)
        println("Core_nonemptiness_check = ", Core_nonemptiness_check)

        # Sh = Shapley_value(length(Players),V_iteration[1:end-1])
        # println()
        # println("Sh = ",round.(Sh, digits=3))

        # # Computing the marginal contributions and the Shapley value in a different way:
        MC_matrix = zeros(2^length(Players),length(Players)) # writing down marginal contributions
        MC_individual = zeros(length(Players)) # write down these specific contributions
        MC_grand_coalition = zeros(length(Players)) # write down these specific contributions
        for pl = 1:length(Players)
            for i = 1:2^length(Players)
                u_with = deepcopy(u[i,:])
                u_without = deepcopy(u_with)
                u_without[pl] = 0 # excluding this player
                row_index_with = findall(all(u .== u_with', dims=2))[1][1]
                row_index_without = findall(all(u .== u_without', dims=2))[1][1]
                MC_matrix[i,pl] = V_iteration[row_index_with] - V_iteration[row_index_without]
                if sum(u_with) == 1 && u_with[pl] ==1 
                    global MC_individual[pl] = MC_matrix[i,pl]
                elseif sum(u_with) == length(Players)
                    global MC_grand_coalition[pl] = MC_matrix[i,pl]
                end
            end
        end

        Sh2 = zeros(length(Players))
        for pl = 1:length(Players)
            for i = 1:2^length(Players)
                # if sum(u[i,:])>=3
                    Sh2[pl] += MC_matrix[i,pl]*(factorial(abs(Int(sum(u[i,:]) - 1))))*factorial(abs(Int(length(Players) - sum(u[i,:]))))/factorial(length(Players))
                # end
            end
        end
        println("Sh2 = ",round.(Sh2, digits=3))
        println("Shapley in % = ",round.(Sh2*100/sum(Sh2), digits=3))
        Sh_softmax, Sh_additive_shift, Sh_discard_negatives = rescale_allocation(Sh2)


        global Sh_iterative_results[Monte_Carlo_iteration,:] = Sh2


        """ Note: V must be negative for `Nucleolus.jl` if allocating profits or savings, and positive if allocating costs. """
        X_Nucleolus = Nucleolus(length(Players),-V_iteration)
        X_Nucleolus = -X_Nucleolus
        println("Nucleolus = ",round.(X_Nucleolus, digits=3))
        println("X_Nucleolus in % = ",round.(X_Nucleolus*100/sum(X_Nucleolus), digits=3))
        X_Nucleolus_softmax, X_Nucleolus_additive_shift, X_Nucleolus_discard_negatives = rescale_allocation(X_Nucleolus)


        # # Computing allocation using the leave-one-out approach (excluding one forecast from the grand coalition):
        global X_leave_one_out = MC_grand_coalition / sum(MC_grand_coalition) * GC_avoided_cost
        println("Leave-one-out = ",X_leave_one_out)
        println("Leave-one-outi in % = ",round.(X_leave_one_out*100/sum(X_leave_one_out), digits=3))
        global X_leave_one_out_softmax, X_leave_one_out_additive_shift, X_leave_one_out_discard_negatives = rescale_allocation(X_leave_one_out)


        if read_predefined_Monte_Carlo_data == false # No Monte Carlo data was given, new generated data will be saved
            iterative_results_for_export[Monte_Carlo_iteration] = Dict(
                :EVA_forecasts => EVA_forecasts,
                :EVA_observations => EVA_observations,

                :Shapley => Sh2,
                :Shapley_softmax => Sh_softmax,
                :Shapley_additive_shift => Sh_additive_shift,
                :Shapley_discard_negatives => Sh_discard_negatives,

                :Nucleolus => X_Nucleolus,
                :Nucleolus_softmax => X_Nucleolus_softmax,
                :Nucleolus_additive_shift => X_Nucleolus_additive_shift,
                :Nucleolus_discard_negatives => X_Nucleolus_discard_negatives,

                :Core_nonemptiness_check => Core_nonemptiness_check,
                :X_leave_one_out => X_leave_one_out,
                :X_leave_one_out_softmax => X_leave_one_out_softmax,
                :X_leave_one_out_additive_shift => X_leave_one_out_additive_shift,
                :X_leave_one_out_discard_negatives => X_leave_one_out_discard_negatives,

                :GC_avoided_cost => V_iteration[length(Players)],
                :GC_DA_dispatch_cost => GC_DA_dispatch_cost,
                :GC_BM_dispatch_cost => GC_BM_dispatch_cost,
                :GC_total_cost => GC_DA_dispatch_cost + GC_BM_dispatch_cost,

                :GC_DA_voltages => GC_DA_voltages,
                :GC_pg_DA => GC_pg_DA,
                :GC_qg_DA => GC_qg_DA,
                :GC_pg_flex_DA_up => GC_pg_flex_DA_up,
                :GC_pg_flex_DA_down => GC_pg_flex_DA_down,
                :GC_qg_flex_DA_up => GC_qg_flex_DA_up,
                :GC_qg_flex_DA_down => GC_qg_flex_DA_down,

                :GC_BM_voltages => GC_BM_voltages,
                :GC_pg_BM_up => GC_pg_BM_up,
                :GC_pg_BM_down => GC_pg_BM_down,
                :GC_qg_BM_up => GC_pg_BM_up,
                :GC_qg_BM_down => GC_pg_BM_down,
                :GC_pg_flex_BM_up => GC_pg_flex_BM_up,
                :GC_pg_flex_BM_down => GC_pg_flex_BM_down,
                :GC_qg_flex_BM_up => GC_qg_flex_BM_up,
                :GC_qg_flex_BM_down => GC_qg_flex_BM_down,

            )
        else
            coalitional_analysis_dict = Dict(
                :coalitional_structure => u,
                :V => V_iteration,
                :MC_matrix => MC_matrix,
                :MC_individual => MC_individual,
                :MC_grand_coalition => MC_grand_coalition,

                :Shapley => Sh2,
                :Shapley_softmax => Sh_softmax,
                :Shapley_additive_shift => Sh_additive_shift,
                :Shapley_discard_negatives => Sh_discard_negatives,


                :Nucleolus => X_Nucleolus,
                :Nucleolus_softmax => X_Nucleolus_softmax,
                :Nucleolus_additive_shift => X_Nucleolus_additive_shift,
                :Nucleolus_discard_negatives => X_Nucleolus_discard_negatives,

                :Core_nonemptiness_check => Core_nonemptiness_check,

                :X_leave_one_out => X_leave_one_out,
                :X_leave_one_out_softmax => X_leave_one_out_softmax,
                :X_leave_one_out_additive_shift => X_leave_one_out_additive_shift,
                :X_leave_one_out_discard_negatives => X_leave_one_out_discard_negatives,
            )
            predefined_Monte_Carlo_data["Monte_Carlo_results_for_export"][Monte_Carlo_iteration][:coalitional_analysis] = coalitional_analysis_dict
        end

        GC.gc() # Force garbage collection?

    else
    println("The coalitional analysis is already performed and saved for this iteration. Skipping it...")
    end

end # <-- end of Monte Carlo iterations

end # elapsed time


EVA_results_errors = EVA_forecasts_iterative_results .- EVA_observations_iterative_results

# # Saving results for all iteration for further analysis:
if read_predefined_Monte_Carlo_data == false # Create a new file with the generated results
    @save "../../results/Monte_Carlo_coalitional_analysis_results_test.jld2" iterative_results_for_export
else # Save to the existing file with Monte Carlo data
    @save "../../results/Monte_Carlo_coalitional_analysis_results_1000_positive_runs.jld2" Monte_Carlo_results_for_export = predefined_Monte_Carlo_data["Monte_Carlo_results_for_export"]
end



""" Plotting of the iterative load analysis starts here """

using Plots
using Plots.PlotMeasures

# fz = 18
fz = 22

plt_Monte_Carlo_coalitional_analysis = plot(
    title = "Monte Carlo coalitional analysis",
    xlabel = "Forecast Error, MW",
    ylabel = "Allocated avoided cost (Shapley value), EUR/h",

    size = (1200,1200), # width and height of the whole plot (in px)
#     aspect_ratio = :equal,

    # xlim = (1,33),
    # ylim = (400, 750),

    xtickfontsize = fz, ytickfontsize = fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,
    legendfontsize = fz,

    # legend = false,
    
    framestyle = :box,
    margin = 10mm,

    minorgrid = :true,
)

for player_i = 1:length(Players)

    scatter!(plt_Monte_Carlo_coalitional_analysis,
            EVA_results_errors[:,player_i],
            Sh_iterative_results[:,player_i],
            label = "EVA "*string(Players[player_i]),
            # color = :black,
            # color = :gray,
            # color = palette(:tab10)[5],
            markersize = 8,
            markerstrokewidth = 0,
            alpha = 0.5
    )
end


display(plt_Monte_Carlo_coalitional_analysis)