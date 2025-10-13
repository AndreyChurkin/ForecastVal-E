"""
This code formulates the OPF optimisation model corresponding to the balancing market dispatch.

This model is not intended to be executed independently.
It is called by higher-level routines that perform the full day-ahead and real-time dispatch optimisation workflow.

Andrey Churkin https://andreychurkin.ru/

"""



if optimize_with == "Gurobi"
    optimizer = Gurobi.Optimizer
    # optimizer = optimizer_with_attributes(() -> Gurobi.Optimizer(), "NonConvex" => 2, "OutputFlag" => 0)
elseif optimize_with == "Ipopt"
    optimizer = Ipopt.Optimizer
else
    println()
    printstyled("Warning: optimizer was not specified correctly!"; color = :magenta)
    println()
end





global OPFmodel_BM = Model(optimizer)

set_silent(OPFmodel_BM) # Suppressing all solver output (including parameter changes messages). Use "unset_silent(model)" to restore output

# set_optimizer_attributes(OPFmodel_BM,"tol" => 1e-05) # Ipopt
set_optimizer_attributes(OPFmodel_BM,"OutputFlag" => 0) # Gurobi
set_optimizer_attributes(OPFmodel_BM,"NonConvex" => 2) # Gurobi
# set_optimizer_attributes(OPFmodel_BM,"MIPGap" => 10^-6) # Gurobi

@variable(OPFmodel_BM, 0 <= pg_BM_up[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["pmax"]) 
@variable(OPFmodel_BM, 0 <= pg_BM_down[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["pmax"])
@variable(OPFmodel_BM, 0 <= qg_BM_up[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["qmax"]) 
@variable(OPFmodel_BM, 0 <= qg_BM_down[i = 1:length(network_data["gen"])] <= network_data["gen"][string(i)]["qmax"]) 

@variable(OPFmodel_BM, p[i = 1:length(network_data["bus"])]) # active load at node i
@variable(OPFmodel_BM, q[i = 1:length(network_data["bus"])]) # reactive load at node i
@variable(OPFmodel_BM, network_data["bus"][string(i)]["vmin"]^2 <= w[i = 1:length(network_data["bus"])] <= network_data["bus"][string(i)]["vmax"]^2) # squared voltage magnitude of node i
@variable(OPFmodel_BM, IL[l = 1:length(network_data["branch"])]) # current magnitute squared
@variable(OPFmodel_BM, pl[l = 1:length(network_data["branch"])]) # active power flow from i to j
@variable(OPFmodel_BM, ql[l = 1:length(network_data["branch"])]) # reactive power flow from i to j

@variable(OPFmodel_BM, 0 <= pg_flex_BM_down[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_p_down"])
@variable(OPFmodel_BM, 0 <= pg_flex_BM_up[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_p_up"])
@variable(OPFmodel_BM, 0 <= qg_flex_BM_down[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_q_down"])
@variable(OPFmodel_BM, 0 <= qg_flex_BM_up[i = 1:length(flex_data["bus"])] <= flex_data["bus"][string(i)]["flex_q_up"])

bus_loads_p = zeros(length(network_data["bus"])) # a vector of nodal actve loads
bus_loads_q = zeros(length(network_data["bus"])) # a vector of nodal reactve loads

for i = 1:length(network_data["load"])
        bus_loads_p[network_data["load"][string(i)]["load_bus"]] = network_data["load"][string(i)]["pd"]
        bus_loads_q[network_data["load"][string(i)]["load_bus"]] = network_data["load"][string(i)]["qd"]
end

@constraint(OPFmodel_BM, [i = 1:length(network_data["bus"])], p[i] == bus_loads_p[i]) # active loads
@constraint(OPFmodel_BM, [i = 1:length(network_data["bus"])], q[i] == bus_loads_q[i]) # reactive loads

# Branch flow equations:
""" Note: DA dispatch of flexible units is included here! """

@constraint(OPFmodel_BM, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                    pl[l]
                    ==
                    p[network_data["branch"][string(l)]["t_bus"]]
                    -
                    pg_flex_BM_up[network_data["branch"][string(l)]["t_bus"]]
                    -
                    DA_dispatch_solution["flex"][string(network_data["branch"][string(l)]["t_bus"])]["pg_flex_up_opt_DA"]
                    +
                    pg_flex_BM_down[network_data["branch"][string(l)]["t_bus"]]
                    +
                    DA_dispatch_solution["flex"][string(network_data["branch"][string(l)]["t_bus"])]["pg_flex_down_opt_DA"]
                    +
                    network_data["branch"][string(l)]["br_r"]*IL[l]
                    +
                    sum(pl[k] for k = 1:1:length(network_data["branch"]) if network_data["branch"][string(k)]["f_bus"] == network_data["branch"][string(l)]["t_bus"])
                    +
                    EVA_data["bus"][string(network_data["branch"][string(l)]["t_bus"])]["load_observation"]
)

@constraint(OPFmodel_BM, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                    ql[l]
                    ==
                    q[network_data["branch"][string(l)]["t_bus"]]
                    -
                    qg_flex_BM_up[network_data["branch"][string(l)]["t_bus"]]
                    +
                    qg_flex_BM_down[network_data["branch"][string(l)]["t_bus"]]
                    +
                    network_data["branch"][string(l)]["br_x"]*IL[l]
                    +
                    sum(ql[k] for k = 1:1:length(network_data["branch"]) if network_data["branch"][string(k)]["f_bus"] == network_data["branch"][string(l)]["t_bus"]) 
)

# Voltage relation:
@constraint(OPFmodel_BM, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                    w[network_data["branch"][string(l)]["t_bus"]]
                    ==
                    w[network_data["branch"][string(l)]["f_bus"]]
                    +
                    (network_data["branch"][string(l)]["br_r"]^2 + network_data["branch"][string(l)]["br_x"]^2)*IL[l]
                    -
                    2*(network_data["branch"][string(l)]["br_r"]*pl[l] + network_data["branch"][string(l)]["br_x"]*ql[l]) 
)

# Power flows and current relation (the exact formulation):
@constraint(OPFmodel_BM, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["br_status"] == 1],
                    pl[l]^2 + ql[l]^2
                    ==
                    IL[l]*w[network_data["branch"][string(l)]["f_bus"]] 
)

# Special conditions for the interface node
""" Note: DA dispatch of the generator is included here! """

if size(gen_line) == ()
        @constraint(OPFmodel_BM, pg_BM_up[1] - pg_BM_down[1] == pl[gen_line] - DA_dispatch_solution["gen"]["p_opt_DA"])
        @constraint(OPFmodel_BM, qg_BM_up[1] - qg_BM_down[1] == ql[gen_line])
else
        @constraint(OPFmodel_BM, pg[1] == sum(pl[line] for line in gen_line) - DA_dispatch_solution["gen"]["p_opt_DA"]) # !? check if this is correct for such cases !?
        @constraint(OPFmodel_BM, qg[1] == sum(ql[line] for line in gen_line) - DA_dispatch_solution["gen"]["p_opt_DA"])
end

if primary_substation_voltage_fixed == true
    @constraint(OPFmodel_BM, w[gen_bus] == 1) # <--- fix the voltage (this may change network's flexibility significantly)
end

if include_line_limits == true # imposing MVA limits for branches
    @constraint(OPFmodel_BM, [l = 1:length(network_data["branch"]); network_data["branch"][string(l)]["rate_a"] >= 10^-4],
                    pl[l]^2 + ql[l]^2 <= network_data["branch"][string(l)]["rate_a"]^2)
end

@objective(OPFmodel_BM, Min, pg_BM_up[1]*TSO_data["cost_BM_up"] + pg_BM_down[1]*TSO_data["cost_BM_down"] 
            + sum(pg_flex_BM_up[i]*flex_data["bus"][string(i)]["cost_BM_up"]
                    + pg_flex_BM_down[i]*flex_data["bus"][string(i)]["cost_BM_down"]
                    for i = 1:length(flex_data["bus"]))
)
