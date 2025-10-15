
"""
This code can be used to test compare different OPF formulations, e.g., models available in PowerModels.jl vs new models

Andrey Churkin https://andreychurkin.ru/

"""


cd(dirname(@__FILE__))
pwd()

using PowerModels, JuMP, Ipopt, Gurobi


# # Select the .m network data file to solve OPF:

# file = "../cases/case33bw (pure).m" # Note: the original case33bw has parsing problems
file = "../cases/case33bw (mod) v1.0.m" # <-- voltage at node 1 fixed to 1.0 pu
# file = "../cases/case33bw (mod) v1.1.m" # <-- case33bw used in the paper "Tracing, Ranking and Valuation of Aggregated DER Flexibility in Active Distribution Networks"




network_data = parse_file(file) # parsing using PowerModels.jl

total_pd = sum(item["pd"] for item in values(network_data["load"]))
total_qd = sum(item["qd"] for item in values(network_data["load"]))

println()
println("Total P load of the network: ",total_pd*network_data["baseMVA"]," MW")
println("Total Q load of the network: ",total_qd*network_data["baseMVA"], "MVAr")
println()



# # Select the optimiser and model formulation (witin PowerModels.jl):

optimizer = Ipopt.Optimizer
# optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter" => 6000)

# optimizer = Gurobi.Optimizer
# optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2)
# optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "NonConvex" => 2, "MIPGap" => 10^-8)

# model = ACPPowerModel # AC power flow Model with polar bus voltage variables
model = ACRPowerModel # AC power flow Model with rectangular bus voltage variables
# model = SOCWRPowerModel # Second-order cone relaxation of bus injection model of AC OPF.



# # Examples of modifying data and model instantiating:

# network_data["load"]["3"]["pd"] = 0.0 # to modify data

network_data["bus"]["1"]["vmax"] = 1.0 # <-- fixed voltage at the primary substation - does not work??
network_data["bus"]["1"]["vmin"] = 1.0
network_data["gen"]["1"]["pmin"] = 0.0 # ??

# pm = instantiate_model(file, model, PowerModels.build_opf)

# print(pm.model)



PowerModels_solution = solve_opf(file, model, optimizer)

vm_solution = zeros(length(PowerModels_solution["solution"]["bus"]))
for bus = 1:length(PowerModels_solution["solution"]["bus"])
      vm_solution[bus] = sqrt(PowerModels_solution["solution"]["bus"][string(bus)]["vr"]^2 + PowerModels_solution["solution"]["bus"][string(bus)]["vi"]^2)
end


PowerModels.print_summary(PowerModels_solution["solution"])



for load_i = 1:length(network_data["load"])
      println("bus #",network_data["load"][string(load_i)]["source_id"][2],", pd = ",network_data["load"][string(load_i)]["pd"],", qd = ",network_data["load"][string(load_i)]["qd"])
end