"""
This code performs the coalitional analysis of EVA forecasts, i.e., forecast valuation.
Marginal contributions of players are estimated for each coalition.

Then, several solution concepts are implemented to allocate the savings among forecast providers, including:
- The Shapley value
- Nucleolus
- Nash bargaining solution
- The Core
- Leave-one-out approach
- Accuracy-adjusted Shapley value

In addition, various rescaled versions of the above allocations can be used to avoid negative payments.

Forecast valuation is performed once in this script, that is, for a single combination of EVA forecasts.

Andrey Churkin https://andreychurkin.ru/

"""


cd(dirname(@__FILE__))
pwd()

using PowerModels, JuMP, Ipopt, Gurobi
using Suppressor


# # Select the solver to use in the OPF formulations:
optimize_with = "Gurobi"
# optimize_with = "Ipopt"

# # Select the .m network data file to solve OPF:

# file = "../cases/case33bw (pure).m" # the original case33bw has parsing problems
# file = "../cases/case33bw (mod) v1.0.m" # <-- voltage at node 1 fixed to 1.0 pu
file = "../cases/case33bw (mod) v1.1.m" # <-- case33bw used in the paper "Tracing, Ranking and Valuation of Aggregated DER Flexibility in Active Distribution Networks"


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
EVA_data["bus"]["22"]["load_observation"] = 0.25/network_data["baseMVA"]
EVA_data["bus"]["22"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

EVA_data["bus"]["32"]["load_forecast"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["32"]["load_observation"] = 0.25/network_data["baseMVA"]
EVA_data["bus"]["32"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

EVA_data["bus"]["14"]["load_forecast"] = 0.2/network_data["baseMVA"]
EVA_data["bus"]["14"]["load_observation"] = 0.25/network_data["baseMVA"]
EVA_data["bus"]["14"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

# EVA_data["bus"]["2"]["load_forecast"] = 0.2/network_data["baseMVA"]
# EVA_data["bus"]["2"]["load_observation"] = 0.25/network_data["baseMVA"]
# EVA_data["bus"]["2"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

# EVA_data["bus"]["3"]["load_forecast"] = 0.2/network_data["baseMVA"]
# EVA_data["bus"]["3"]["load_observation"] = 0.25/network_data["baseMVA"]
# EVA_data["bus"]["3"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

# EVA_data["bus"]["4"]["load_forecast"] = 0.2/network_data["baseMVA"]
# EVA_data["bus"]["4"]["load_observation"] = 0.25/network_data["baseMVA"]
# EVA_data["bus"]["4"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

# EVA_data["bus"]["5"]["load_forecast"] = 0.2/network_data["baseMVA"]
# EVA_data["bus"]["5"]["load_observation"] = 0.25/network_data["baseMVA"]
# EVA_data["bus"]["5"]["status_quo_forecast"] = 0.1/network_data["baseMVA"]

"""
Note: even if EVAs are not included as players, they will still impact the DA and BM dispatches, affecting the value of data from other players
If not including an EVA in the coalitional game, delete its data from EVA_data

"""

# # Define the players of a cooperative game: names of EVAs (their buses), and their order:
# # One player:
# Players = [14]

# # 2 players:
# Players = [32 14]
# Players = [22 32]

# # Three players:
Players = [22 32 14]
# Players = [22 32 5]
# Players = [22 2 3]
# Players = [2 3 4]


# # 4 players:
# Players = [22 32 14 5]
# Players = [2 3 4 5]



# # ------------ Data valuation starts here ------------

# # Functions for defining the coalitional structure and calculating the Shapley value:
include("Shapley_matrix.jl")
include("Shapley_value.jl")
include("allocation_rescaling.jl")
# # Function to check whether the Core of a cooperative game is non-empty (returns true if non-empty, false otherwise):
include("Core_nonemptiness_check.jl")
# # Function to calculate the Nucleolus:
include("Nucleolus.jl")
# # Function to calculate the Nash bargaining solution:
include("Nash_Bargaining.jl")
# # Function to calculate the accuracy-adjusted Shapley value:
include("Accuracy_adjusted_Shapley.jl")

u = Shapley_matrix(length(Players)) # get the coalitional structure
V = zeros(2^length(Players)) # the characteristic function

u_buses = zeros(size(u,1), N) # relation between buses and coalitions: EVAs at which buses belong to which coalitions

for coalition = 1:2^length(Players)
      for player = 1:length(Players)
            if u[coalition,player] == 1
                  global u_buses[coalition,Players[player]] = 1
            end
      end
end



# # First, we need to estimate the cost with no EVAs forecasts (only with DSO's status quo forecast):

println()
println("Solving the day-ahead and balancing dispatches for the case with no EVA forecasts")

@suppress begin
      include("DA_dispatch_without_forecasts.jl") # <-- day-ahead dispatch formulation
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
      include("BM_dispatch.jl") # <-- balancing dispatch formulation
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

println("Total cost (DA + BM) = ",total_DA_and_BM_cost_without_forecast, " EUR/h")
println()


for cl = 1:2^length(Players)
      global coalition = cl
      println()
      println("---------- Coalition # ",coalition, " ----------")
      println("Players in coalition = ",u[coalition,:])
      println()

      # # First, the day-ahead dispatch is formulated and solved:

      @suppress begin
            include("DA_dispatch_for_coalition.jl") # <-- day-ahead dispatch formulation
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

      println("Primary substation DA supply = ",round(value.(pg[gen_bus])*network_data["baseMVA"],digits=3)," MW, ",round(value.(qg[gen_bus])*network_data["baseMVA"],digits=3)," MVAr")
      println("Flexibility DA dispatch: upward = ",round(value.(pg_flex_up)[7]*network_data["baseMVA"],digits=3)," MW, downward = ",round(value.(pg_flex_down)[7]*network_data["baseMVA"],digits=3)," MVAr")
      println("DA dispatch cost = ",round(objective_value(OPFmodel_DA)*network_data["baseMVA"],digits=3)," EUR/h")
      println()

      # # Then, the balancing market dispatch is formulated and solved:

      @suppress begin
            include("BM_dispatch.jl") # <-- balancing dispatch formulation
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

      println("Primary substation balancing:")
      println("P_up = ",round(value.(pg_BM_up[gen_bus])*network_data["baseMVA"],digits=3)," MW, Q_up = ",round(value.(qg_BM_up[gen_bus])*network_data["baseMVA"],digits=3)," MVAr")
      println("P_down = ",round(value.(pg_BM_down[gen_bus])*network_data["baseMVA"],digits=3)," MW, Q_down = ",round(value.(qg_BM_down[gen_bus])*network_data["baseMVA"],digits=3)," MVAr")
      println("Flexibility balancing dispatch: upward = ",round(value.(pg_flex_BM_up)[7]*network_data["baseMVA"],digits=3)," MW, downward = ",round(value.(pg_flex_BM_down)[7]*network_data["baseMVA"],digits=3)," MVAr")
      println("Balancing dispatch cost = ",round(objective_value(OPFmodel_BM)*network_data["baseMVA"],digits=3)," EUR/h")

      min_voltage_BM, min_voltage_idx_BM = findmin(sqrt.(value.(w)))

      println("Lowest voltage observed = ",round(min_voltage_BM,digits=3)," kV, at bus #",min_voltage_idx_BM)

      println()

      println("Total cost (DA + BM) = ",round(DA_dispatch_solution["total_DA_cost"]["EUR/h"] + BM_dispatch_solution["total_BM_cost"]["EUR/h"],digits=3), " EUR/h")
      
      global coalition_avoided_cost = total_DA_and_BM_cost_without_forecast - (DA_dispatch_solution["total_DA_cost"]["EUR/h"] + BM_dispatch_solution["total_BM_cost"]["EUR/h"])
      
      println("Avoided cost due to the coalition's aggregate forecast = ", round(coalition_avoided_cost,digits=3)," EUR/h")
      println()
      
      global V[coalition] = coalition_avoided_cost

end

println("---------- Coalitional analysis completed ----------")

GC_avoided_cost = V[length(Players)] # avoided cost of the grand coalition

println()
println("Avoided cost of the Grand Coalition = ", round(GC_avoided_cost,digits=3))


global Core_nonemptiness_check = Core_nonemptiness(n = length(Players), coalitional_structure = u, V = V)



Sh = Shapley_value(length(Players),V[1:end-1])
println()
println("Sh = ",round.(Sh, digits=3))

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
        MC_matrix[i,pl] = V[row_index_with] - V[row_index_without]
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

# # Rescaling the Shapley value:
Sh_softmax, Sh_additive_shift, Sh_discard_negatives = rescale_allocation(Sh)
println("Sh_softmax = ",round.(Sh_softmax, digits=3))
println("Sh_softmax in % = ",round.(Sh_softmax*100/sum(Sh_softmax), digits=3))
println("Sh_additive_shift = ",round.(Sh_additive_shift, digits=3))
println("Sh_additive_shift in % = ",round.(Sh_additive_shift*100/sum(Sh_additive_shift), digits=3))
println("Sh_discard_negatives = ",round.(Sh_discard_negatives, digits=3))
println("Sh_discard_negatives in % = ",round.(Sh_discard_negatives*100/sum(Sh_discard_negatives), digits=3))

""" Note: V must be negative for `Nucleolus.jl` if allocating profits or savings, and positive if allocating costs. """
X_Nucleolus = Nucleolus(length(Players),-V)
X_Nucleolus = -X_Nucleolus
println("X_Nucleolus = ",round.(X_Nucleolus, digits=3))
println("X_Nucleolus in % = ",round.(X_Nucleolus*100/sum(X_Nucleolus), digits=3))



# # # Option 1: Defining the disagreement points as zeros (player gets nothing if not sharing forecast)
# global d = zeros(length(Players))

# # Option 1: Defining the disagreement points at the contributions to the grand coalition v(N)-v(N∖{i})
global d = MC_grand_coalition

""" Note: V is expected to be positive for `Nash_Bargaining.jl` when allocating profits or savings """
# X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(length(Players), d; V = V[1:end-1], objective_formulation = "log")
global X_Nash_bargaining, Nash_bargaining_model = Nash_bargaining(length(Players), d; total_value = V[length(Players)], objective_formulation = "log")

println("X_Nash_bargaining = ",round.(X_Nash_bargaining, digits=3))
println("X_Nash_bargaining in % = ",round.(X_Nash_bargaining*100/sum(X_Nash_bargaining), digits=3))



# # Computing allocation using the leave-one-out approach (excluding one forecast from the grand coalition):
X_leave_one_out = MC_grand_coalition / sum(MC_grand_coalition) * GC_avoided_cost
X_leave_one_out_softmax, X_leave_one_out_additive_shift, X_leave_one_out_discard_negatives = rescale_allocation(X_leave_one_out)
println()
println("X_leave_one_out = ",round.(X_leave_one_out, digits=3))
println("X_leave_one_out in % = ",round.(X_leave_one_out*100/sum(X_leave_one_out), digits=3))
println("X_leave_one_out_softmax = ",round.(X_leave_one_out_softmax, digits=3))
println("X_leave_one_out_softmax in % = ",round.(X_leave_one_out_softmax*100/sum(X_leave_one_out_softmax), digits=3))
println("X_leave_one_out_additive_shift = ",round.(X_leave_one_out_additive_shift, digits=3))
println("X_leave_one_out_additive_shift in % = ",round.(X_leave_one_out_additive_shift*100/sum(X_leave_one_out_additive_shift), digits=3))
println("X_leave_one_out_discard_negatives = ",round.(X_leave_one_out_discard_negatives, digits=3))
println("X_leave_one_out_discard_negatives in % = ",round.(X_leave_one_out_discard_negatives*100/sum(X_leave_one_out_discard_negatives), digits=3))



# # Calculating the accuracy-adjusted Shapley value:
EVA_forecast_errors = zeros(length(Players))
for player = 1:length(Players)
      forecast = EVA_data["bus"][string(Players[player])]["load_forecast"]
      observation = EVA_data["bus"][string(Players[player])]["load_observation"]
      error = forecast - observation
      EVA_forecast_errors[player] = error
end

accuracy_adjusted_Shapley, error_scores, only_accuracy_based_payments = AA_Shapley(
      V_GC = GC_avoided_cost,
      original_Shapley = Sh2,
      errors = EVA_forecast_errors,
      α = 0.3, # fraction of value allocated based on accuracy
      ρ = 2.0,
      τ = 0.05
)

println()
println("accuracy_adjusted_Shapley = ", round.(accuracy_adjusted_Shapley, digits=3))
println("accuracy_adjusted_Shapley in % = ", round.(accuracy_adjusted_Shapley*100/sum(accuracy_adjusted_Shapley), digits=3))
println("error_scores = ", round.(error_scores, digits=3))
println("only_accuracy_based_payments = ", round.(only_accuracy_based_payments, digits=3))
println("only_accuracy_based_payments in % = ", round.(only_accuracy_based_payments*100/sum(only_accuracy_based_payments), digits=3))




""" Plotting of the marginal contributions and the Shapley value starts here """


using Plots, Plots.PlotMeasures
using StatsPlots

# fz = 20
fz = 22

# Let’s not plot the zero values for coalitions where player is not part of:
MC_matrix_no_zeros = zeros(Int((2^length(Players))/2), length(Players))

for pl = 1:Int(length(Players))
MC_matrix_no_zeros[:,pl] = collect(x for x in MC_matrix[:,pl] if x != 0)
end


plt = violin(collect(1:length(Players))',
      # MC_matrix,
      MC_matrix_no_zeros,

    #   alpha=0.25,
      alpha=0.2,
    #   alpha=1,
      # alpha=0.75,
      lw = 3,
      # linecolor =
      # color=palette(:tab10)[1],
      # color=palette(:tab10)[2],
      # color=palette(:tab10)[3],
      color = :grey,
    #   color=palette(:tab10)[6],
      # color = RGB(flex_color[1],flex_color[2],flex_color[3]),
      # linealpha = 0.7,
      # linealpha = 0.0,
      # linewidth=0.0,
      fontfamily = "Courier",
      size = (1200,1200),
    #   xlim=(0,nPl+1),
      # ylim = (0,maximum(time_records_loops)),
    #   ylim = (0,160),
      xlabel = "EVAs (players)", 
      ylabel = "Contribution to avoided cost, EUR/h",
      xtickfontsize = fz, ytickfontsize = fz, xguidefontsize = fz, yguidefontsize = fz, legendfontsize = fz,
      foreground_color_legend = nothing, 
      legend = false,
      framestyle = :box, margin = 20mm, left_margin = 50mm, yminorgrid = :true,
      xticks = (1:length(Players), string.(1:length(Players))),
      # aspect_ratio=:equal
      # fill_z = 1 - (f)/F,
      # color = palette([RGB(0,0,0), RGB(1,1,1)], F)
)

for pl = 1:length(Players)

    scatter!(plt,
            # ones(size(MC_matrix)[1])*pl, MC_matrix[:,pl],
            ones(size(MC_matrix_no_zeros)[1])*pl, MC_matrix_no_zeros[:,pl],

            # markersize = 20,
            markersize = 15,
            # markershape = :cross,
            # markershape = :hline,
            markershape = :circle,
            # markercolor = :grey
            markercolor = RGB(0.0,0.0,0.0),
            markerstrokewidth = 0,
            alpha = 0.5
    )

#     scatter!(plt,[pl], [MC_individual[pl]],
#     markersize = 15,
#     # markersize = 5,
#     markershape = :cross,
#     # markershape = :hline,
#     # markershape = :circle,
#     # markercolor = :red,
#     markerstrokewidth = 2,
#     markercolor = RGB(227/255,30/255,36/255)
#     )

#     scatter!(plt,[pl], [MC_grand_coalition[pl]],
#     markersize = 15,
#     # markersize = 5,
#     markershape = :x,
#     markerstrokewidth = 2,
#     # markershape = :hline,
#     # markershape = :circle,
#     # markercolor = :red
#     markercolor = RGB(227/255,30/255,36/255)
#     )

    scatter!(plt,[pl], [Sh2[pl]],
    markersize = 15,
    # markersize = 5,
    # markershape = :xcross,
    # markershape = :hline,
    markershape = :diamond,
    # markercolor = :red,
    # markerstrokecolor = :red,
    # markerstrokewidth = 1.2,
    markercolor = RGB(227/255,30/255,36/255),
    label = "Shapley value"
    )
end

display(plt)

# savefig("../results/violin_plot.svg")
# savefig("../results/violin_plot.png")
# savefig("../results/violin_plot.pdf")