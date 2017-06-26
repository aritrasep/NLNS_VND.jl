#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

#####################################################################
# Perturbations based on Nodes                                      #
#####################################################################

#####################################################################
# Perturbation 1                                                    #
#####################################################################

@inbounds function Perturbation_1(data::BikeShareData,solution::BikeShareSolution,repairing_heuristic::Symbol)
	node_rank_list,max_perturbations = sorted_node_rank_list(data,solution)
	perturbations::Int64 = 1
	while perturbations <= (max_perturbations[1]*15)/data.N
		if perturbations%2 != 0
			ind = node_rank_list[1][1:round(Int64,(perturbations)/2)]
		else
			ind = node_rank_list[1][1:max_perturbations[1]+1-round(Int64,(perturbations)/2)]
		end
		if length(ind) > 1
			sort!(ind)
		end
		tmp_solution = deepcopy(solution)
		deleteat!(tmp_solution.Tour[1],ind)
		(tmp_solution.Tour[1],tmp_solution.Instructions[1]) = TOUR_REPAIR(data,tmp_solution.Tour[1],repairing_heuristic)
		VND(data,tmp_solution)
		if tmp_solution.Makespan < solution.Makespan
			solution = deepcopy(tmp_solution)
			node_rank_list,max_perturbations = sorted_node_rank_list(data,solution)
			perturbations = 1
		else
			perturbations += 1
		end
	end
	solution
end

#####################################################################
# Perturbation 2                                                    #
#####################################################################

@inbounds function Perturbation_2(data::BikeShareData,solution::BikeShareSolution,repairing_heuristic::Symbol)
	best_solution::BikeShareSolution = deepcopy(solution)
	node_rank_list,max_perturbations = sorted_node_rank_list(data,best_solution)
	ind = node_rank_list[1][1:max_perturbations[1]]
	i = 2
	while i <= length(ind)
		if ind[i] > ind[i-1]
			deleteat!(ind,i)
		else 
			i += 1
		end 
	end
	i = 1
	while i <= length(ind) 
		tmp_solution = deepcopy(best_solution)
		tmp_solution.Tour[1] = [tmp_solution.Tour[1][1:ind[i]];1]
		(tmp_solution.Tour[1],tmp_solution.Instructions[1]) = TOUR_REPAIR(data,tmp_solution.Tour[1],repairing_heuristic)
		VND(data,tmp_solution)
		if tmp_solution.Makespan < best_solution.Makespan
			best_solution = deepcopy(tmp_solution)
			node_rank_list,max_perturbations = sorted_node_rank_list(data,best_solution)
			ind = node_rank_list[1][1:max_perturbations[1]]
			i = 2
			while i <= length(ind)
				if ind[i] > ind[i-1]
					deleteat!(ind,i)
				else 
					i += 1
				end 
			end
			i = 1
		else
			i += 1
		end
	end
	best_solution
end

#####################################################################
# Perturbation 3                                                    #
#####################################################################

@inbounds function Perturbation_3(data::BikeShareData,solution::BikeShareSolution,repairing_heuristic::Symbol)
	best_solution::BikeShareSolution = deepcopy(solution)
	node_rank_list,max_perturbations = sorted_node_rank_list(data,best_solution)
	ind = node_rank_list[1][1:max_perturbations[1]]
	i = 2
	while i <= length(ind)
		if ind[i] < ind[i-1]
			deleteat!(ind,i)
		else 
			i += 1
		end 
	end
	i = 1
	while i <= length(ind) 
		tmp_solution = deepcopy(best_solution)
		tmp_solution.Tour[1] = [1;tmp_solution.Tour[1][ind[i]:end]]
		(tmp_solution.Tour[1],tmp_solution.Instructions[1]) = TOUR_REPAIR(data,tmp_solution.Tour[1],repairing_heuristic)
		VND(data,tmp_solution)
		if tmp_solution.Makespan < best_solution.Makespan
			best_solution = deepcopy(tmp_solution)
			node_rank_list,max_perturbations = sorted_node_rank_list(data,best_solution)
			ind = node_rank_list[1][1:max_perturbations[1]]
			i = 2
			while i <= length(ind)
				if ind[i] < ind[i-1]
					deleteat!(ind,i)
				else 
					i += 1
				end 
			end
			i = 1
		else
			i += 1
		end
	end
	best_solution
end

#####################################################################
# Perturbation 4                                                    #
#####################################################################

@inbounds function Perturbation_4(data::BikeShareData,solution::BikeShareSolution,repairing_heuristic::Symbol)
	node_rank_list,max_perturbations = sorted_node_rank_list(data,solution)
	node_rank_list[1] = reverse(node_rank_list[1])
	max_perturbations[1] = length(node_rank_list[1]) - max_perturbations[1]
	perturbations::Int64 = 1
	while perturbations <= (max_perturbations[1]*15)/data.N
		if perturbations%2 != 0
			ind = node_rank_list[1][1:round(Int64,(perturbations)/2)]
		else
			ind = node_rank_list[1][1:max_perturbations[1]+1-round(Int64,div(perturbations,2))]
		end
		if length(ind) > 1
			sort!(ind)
		end
		tmp_solution = deepcopy(solution)
		deleteat!(tmp_solution.Tour[1],ind)
		(tmp_solution.Tour[1],tmp_solution.Instructions[1]) = TOUR_REPAIR(data,tmp_solution.Tour[1],repairing_heuristic)
		VND(data,tmp_solution)
		if tmp_solution.Makespan < solution.Makespan
			solution = deepcopy(tmp_solution)
			node_rank_list,max_perturbations = sorted_node_rank_list(data,solution)
			node_rank_list[1] = reverse(node_rank_list[1])
			max_perturbations[1] = length(node_rank_list[1]) - max_perturbations[1]
			perturbations = 1
		else
			perturbations += 1
		end
	end
	solution
end

#####################################################################
# Deep Exploration of Neighborhoods for Single Vehicle              #
#####################################################################

@inbounds function DEEP_EXPLORATION(data::BikeShareData,solution::BikeShareSolution)
	best_solution::BikeShareSolution = deepcopy(solution)
	perturbation_heuristics::Vector{Symbol}, repairing_heuristics::Vector{Symbol} = [:Perturbation_1,:Perturbation_4], [:Nearest_Neighbor_1,:Nearest_Neighbor_2,:Nearest_Neighbor_3]
	shuffle!(perturbation_heuristics); shuffle!(repairing_heuristics)
	improvement_list::Vector{Float64} = [best_solution.Makespan]
	exit = false
	while exit == false
		for perturbation_heuristic in perturbation_heuristics, repairing_heuristic in repairing_heuristics
			tmp_best_solution = eval(perturbation_heuristic)(data,best_solution,repairing_heuristic)
			if best_solution.Makespan > tmp_best_solution.Makespan
				best_solution = deepcopy(tmp_best_solution)
			end
			push!(improvement_list,best_solution.Makespan)
			if length(improvement_list) > (length(perturbation_heuristics)*length(repairing_heuristics)) && improvement_list[end] == improvement_list[end-(length(perturbation_heuristics)*length(repairing_heuristics))]
				exit = true
				break
			end
		end
		if improvement_list[end] != solution.Makespan
			perturbation_heuristics = [:Perturbation_2, :Perturbation_3, :Perturbation_1, :Perturbation_4]
		end
	end
	best_solution
end
