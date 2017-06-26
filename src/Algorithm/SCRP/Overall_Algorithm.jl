#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

include("Construction_and_Repairing_Operators.jl")
include("Local_Search_Operators1.jl")
include("Local_Search_Operators2.jl")
include("Variable_Neighborhood_Descent.jl")
include("Perturbation_Operators.jl")

#####################################################################
# NLNS+VND for Single Vehicle                                       #
#####################################################################

@inbounds function nlns_vnd_1V(data::BikeShareData)
	t0::Float64 = time()
	solution::BikeShareSolution = INITIAL_SOLUTION_CREATION(data)
	solution = nlns_vnd_1V(data,solution)
	t1::Float64 = time()
	solution.Computing_Time = round(t1-t0,6)
	solution
end

@inbounds function nlns_vnd_1V(data::BikeShareData,solution::BikeShareSolution)
	DEEP_EXPLORATION(data,solution)
end

#####################################################################
# NLNS+VND for Multiple Vehicles                                    #
#####################################################################

@inbounds function nlns_vnd_MV(data::BikeShareData)
	t0::Float64 = time()
	solution::BikeShareSolution = INITIAL_SOLUTION_CREATION(data)
	solution = nlns_vnd_MV(data,solution)
	t1::Float64 = time()
	solution.Computing_Time = round(t1-t0,6)
	solution
end

@inbounds function nlns_vnd_MV(data::BikeShareData,solution::BikeShareSolution)
	improvement_list::Vector{Float64} = [solution.Makespan]
	split_data::Vector{BikeShareData}, split_solution::Vector{BikeShareSolution} = BikeShareData[], BikeShareSolution[]
	while true
		split_data,split_solution = split_bikeshare_data_and_solution(data,solution)
		for i in 1:data.V
			split_solution[i] = nlns_vnd_1V(split_data[i],split_solution[i])
		end
		solution = combine_split_bikeshare_solution(split_solution)
		VND(data,solution)
		push!(improvement_list,solution.Makespan)
		if improvement_list[end] == improvement_list[end-1]
			break
		end
	end
	solution
end

@inbounds function split_bikeshare_data_and_solution(data::BikeShareData,solution::BikeShareSolution)
	split_data::Vector{BikeShareData}, split_solution::Vector{BikeShareSolution}, Tour::Vector{Vector{Int64}}, Instructions::Vector{Vector{Int64}} = BikeShareData[], BikeShareSolution[], Vector{Int64}[], Vector{Int64}[]
	for i in 1:data.V
		push!(split_data,deepcopy(data))
		split_data[i].O, split_data[i].V = zeros(Int64,data.N), 1		
		for j in 1:length(solution.Tour[i])
			split_data[i].O[solution.Tour[i][j]] += solution.Instructions[i][j]
		end
		Tour, Instructions = Vector{Int64}[], Vector{Int64}[]
		push!(Tour,solution.Tour[i]); push!(Instructions, solution.Instructions[i])
		push!(split_solution,BikeShareSolution(solution.Rebalancing_Times[i],[solution.Rebalancing_Times[i]],Tour,Instructions,solution.Rank_List,0.0))
	end
	split_data,split_solution
end

@inbounds function combine_split_bikeshare_solution(split_solution::Vector{BikeShareSolution})
	Rebalancing_Times::Vector{Float64}, Tour::Vector{Vector{Int64}}, Instructions::Vector{Vector{Int64}} = Float64[], Vector{Int64}[], Vector{Int64}[]
	for i in 1:length(split_solution)
		push!(Rebalancing_Times,split_solution[i].Makespan); push!(Tour,split_solution[i].Tour[1]); push!(Instructions,split_solution[i].Instructions[1])
	end
	BikeShareSolution(maximum(Rebalancing_Times),Rebalancing_Times,Tour,Instructions,split_solution[1].Rank_List,0.0)
end

#####################################################################
# Overall Algorithm                                                 #
#####################################################################

#####################################################################
# NLNS+VND                                                          #
#####################################################################

@inbounds function nlns_vnd(data::BikeShareData)
	@match data.V begin
		1 => nlns_vnd_1V(data)
		_ => nlns_vnd_MV(data)
	end
end
