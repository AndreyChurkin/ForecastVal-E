
"""
This code can be used to analyse a distribution network case with a DistFlow OPF formulation.

Andrey Churkin https://andreychurkin.ru/

"""


cd(dirname(@__FILE__))
pwd()

using PowerModels, JuMP, Ipopt, Gurobi


# # Select the .m network data file to solve OPF:

# file = "../cases/case33bw (pure).m" # Note: the original case33bw has parsing problems
# file = "../cases/case33bw (mod) v1.0.m" # <-- voltage at node 1 fixed to 1.0 pu
file = "../cases/case33bw (mod) v1.1.m" # <-- case33bw used in the paper "Tracing, Ranking and Valuation of Aggregated DER Flexibility in Active Distribution Networks"




network_data = parse_file(file)

L = length(network_data["branch"])
N = length(network_data["bus"])

total_pd = sum(item["pd"] for item in values(network_data["load"]))
total_qd = sum(item["qd"] for item in values(network_data["load"]))

println()
println("Total P load of the network: ",total_pd*network_data["baseMVA"]," MW")
println("Total Q load of the network: ",total_qd*network_data["baseMVA"], "MVAr")
println()

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



# # Use this part to introduce flexible generators:

# creating a dictionary for flexible units' data:
flex_data = Dict("bus" =>
            Dict(string(i) =>
                    Dict("bus_i" => i, "flex_p_up" => 0.0, "flex_p_down" => 0.0, "flex_q_up" => 0.0, "flex_q_down" => 0.0, "r" => 1.0, "cost" => 0.0)
                    for i = 1:N ) 
)

# # # defining flexible units in the network:
# flex_data["bus"]["17"]["flex_p_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["17"]["flex_p_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["17"]["flex_q_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["17"]["flex_q_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["17"]["r"] = 0.92 # reliability level -- not needed for now
# flex_data["bus"]["17"]["cost"] = 350 # activation cost $/MWh

# flex_data["bus"]["36"]["flex_p_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["36"]["flex_p_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["36"]["flex_q_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["36"]["flex_q_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["36"]["r"] = 0.92 # reliability level -- not needed for now
# flex_data["bus"]["36"]["cost"] = 300 # activation cost $/MWh

# # # Adding more units:
# flex_data["bus"]["12"]["flex_p_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["12"]["flex_p_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["12"]["flex_q_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["12"]["flex_q_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["12"]["r"] = 0.92 # reliability level -- not needed for now
# flex_data["bus"]["12"]["cost"] = 325 # activation cost $/MWh

# flex_data["bus"]["24"]["flex_p_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["24"]["flex_p_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["24"]["flex_q_up"] = 1/network_data["baseMVA"]
# flex_data["bus"]["24"]["flex_q_down"] = 1/network_data["baseMVA"]
# flex_data["bus"]["24"]["r"] = 0.92 # reliability level -- not needed for now
# flex_data["bus"]["24"]["cost"] = 375 # activation cost $/MWh



# # Defining the DistFlow OPF formulation with upward and downward regulations

function build_flex_DistFlow(network_data, flex_data) # Note: OPF formulation with upward and downward regulations

    # optimizer = Ipopt.Optimizer
    optimizer = Gurobi.Optimizer

    global OPFmodel = Model(optimizer)

    # set_optimizer_attributes(OPFmodel,"tol" => 1e-05) # Ipopt
    set_optimizer_attributes(OPFmodel,"NonConvex" => 2) # Gurobi
    set_optimizer_attributes(OPFmodel,"OutputFlag" => 1) # Gurobi
    # set_optimizer_attributes(OPFmodel,"MIPGap" => 10^-6) # Gurobi

    @variable(OPFmodel, network_data["gen"][string(i)]["pmin"] <= pg[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["pmax"]) # generation power (real) at node i
    @variable(OPFmodel, network_data["gen"][string(i)]["qmin"] <= qg[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["qmax"]) # generation power (imaginary) at node i
    @variable(OPFmodel, p[i = 1:length(network_data["bus"])]) # active load at node i
    @variable(OPFmodel, q[i = 1:length(network_data["bus"])]) # reactive load at node i
    @variable(OPFmodel, network_data["bus"][string(i)]["vmin"]^2 <= w[i = 1:length(network_data["bus"])] <= network_data["bus"][string(i)]["vmax"]^2) # squared voltage magnitude of node i
    @variable(OPFmodel, IL[l = 1:length(network_data["branch"])]) # current magnitute squared
    @variable(OPFmodel, pl[l = 1:length(network_data["branch"])]) # active power flow from i to j
    @variable(OPFmodel, ql[l = 1:length(network_data["branch"])]) # reactive power flow from i to j

    @variable(OPFmodel, 0 <= pg_flex_down[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_p_down"])
    @variable(OPFmodel, 0 <= pg_flex_up[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_p_up"])
    @variable(OPFmodel, 0 <= qg_flex_down[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_q_down"])
    @variable(OPFmodel, 0 <= qg_flex_up[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_q_up"])

    global pg, qg, p, q, w, IL, pl, ql, pg_flex_up, pg_flex_down, qg_flex_up, qg_flex_down # <-- make variables accessible outside of the function

    bus_loads_p = zeros(length(network_data["bus"])) # a vector of nodal actve loads
    bus_loads_q = zeros(length(network_data["bus"])) # a vector of nodal reactve loads

    for i = 1:length(network_data["load"])
          bus_loads_p[network_data["load"][string(i)]["load_bus"]] = network_data["load"][string(i)]["pd"]
          bus_loads_q[network_data["load"][string(i)]["load_bus"]] = network_data["load"][string(i)]["qd"]
    end

    @constraint(OPFmodel, [i = 1:length(network_data["bus"])], p[i] == bus_loads_p[i]) # active loads
    @constraint(OPFmodel, [i = 1:length(network_data["bus"])], q[i] == bus_loads_q[i]) # reactive loads

    # Branch flow equations:
    @constraint(OPFmodel, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                      pl[l]
                      ==
                      p[network_data["branch"][string(l)]["t_bus"]]
                      -
                      pg_flex_up[network_data["branch"][string(l)]["t_bus"]]
                      +
                      pg_flex_down[network_data["branch"][string(l)]["t_bus"]]
                      +
                      network_data["branch"][string(l)]["br_r"]*IL[l]
                      +
                      sum(pl[k] for k = 1:1:length(network_data["branch"]) if network_data["branch"][string(k)]["f_bus"] == network_data["branch"][string(l)]["t_bus"]) 
    )

    @constraint(OPFmodel, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                      ql[l]
                      ==
                      q[network_data["branch"][string(l)]["t_bus"]]
                      -
                      qg_flex_up[network_data["branch"][string(l)]["t_bus"]]
                      +
                      qg_flex_down[network_data["branch"][string(l)]["t_bus"]]
                      +
                      network_data["branch"][string(l)]["br_x"]*IL[l]
                      +
                      sum(ql[k] for k = 1:1:length(network_data["branch"]) if network_data["branch"][string(k)]["f_bus"] == network_data["branch"][string(l)]["t_bus"]) 
    )

    # Voltage relation:
    @constraint(OPFmodel, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                      w[network_data["branch"][string(l)]["t_bus"]]
                      ==
                      w[network_data["branch"][string(l)]["f_bus"]]
                      +
                      (network_data["branch"][string(l)]["br_r"]^2 + network_data["branch"][string(l)]["br_x"]^2)*IL[l]
                      -
                      2*(network_data["branch"][string(l)]["br_r"]*pl[l] + network_data["branch"][string(l)]["br_x"]*ql[l]) 
    )

    # Power flows and current relation (the exact formulation):
    @constraint(OPFmodel, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                      pl[l]^2 + ql[l]^2
                      ==
                      IL[l]*w[network_data["branch"][string(l)]["f_bus"]] 
    )

    # Special conditions for the interface node
    if size(gen_line) == ()
          @constraint(OPFmodel, pg[1] == pl[gen_line])
          @constraint(OPFmodel, qg[1] == ql[gen_line])
    else
          @constraint(OPFmodel, pg[1] == sum(pl[line] for line in gen_line)) # !? check if this is correct for such cases !?
          @constraint(OPFmodel, qg[1] == sum(ql[line] for line in gen_line))
    end

    if primary_substation_voltage_fixed == true
        @constraint(OPFmodel, w[gen_bus] == 1) # <--- fix the voltage (this may change network's flexibility significantly)
    end

    if include_line_limits == true # imposing MVA limits for branches
        @constraint(OPFmodel, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["rate_a"] >= 10^-4],
                        pl[l]^2 + ql[l]^2 <= network_data["branch"][string(l)]["rate_a"]^2)
    end
end



# # Solving the OPF:

build_flex_DistFlow(network_data, flex_data)

@constraint(OPFmodel, [i = 1:length(network_data["bus"])], pg_flex_up[i] == 0) # switch off all units
@constraint(OPFmodel, [i = 1:length(network_data["bus"])], pg_flex_down[i] == 0)
@constraint(OPFmodel, [i = 1:length(network_data["bus"])], qg_flex_up[i] == 0)
@constraint(OPFmodel, [i = 1:length(network_data["bus"])], qg_flex_down[i] == 0)

@objective(OPFmodel, Min, pg[1])
optimize!(OPFmodel)

pl_results_0 = value.(pl)
ql_results_0 = value.(ql)
w_results_0 = value.(w);  v_results_0 = sqrt.(w_results_0)
pg_results_0 = value.(pg[1])
qg_results_0 = value.(qg[1])



# # Plotting the distribution of voltages:

using Plots
using Plots.PlotMeasures

fz = 18
# fz = 22

plt_v_profile = plot(
    title = "Voltage profile",
    xlabel = "bus #",
    ylabel = "Voltage magnitude, p.u.",

#     size = (2000,1000), # width and height of the whole plot (in px)
      size = (2000,500), # width and height of the whole plot (in px)

    # aspect_ratio = :equal,

    xlim = (1,33),
    ylim = (0.85, 1.15),

    xtickfontsize = fz, ytickfontsize = fz,
    fontfamily = "Courier", 
    titlefontsize = fz,
    xguidefontsize = fz,
    yguidefontsize = fz,
    legendfontsize = fz,

    legend = false,
    
    framestyle = :box,
    margin = 20mm,

    minorgrid = :true,
)

plot!(plt_v_profile,
    collect(1:N),
    v_results_0,
    label = "Voltage magnitude, p.u.",
    color = :black,
    w = 8
)

plot!(plt_v_profile,
    [0, 35], [0.9, 0.9], fillrange=[1.1, 1.1], 
    fillalpha = 0.15, 
    fillcolor = palette(:tab10)[3], 
    label="Voltage limit", 
    linewidth = 0
    )

display(plt_v_profile)

savefig("../results/plt_v_profile.png")
# savefig("../results/plt_v_profile.pdf")
# savefig("../results/plt_v_profile.svg")