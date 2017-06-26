#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

#####################################################################
# Functions for computing Rank List                                 #
#####################################################################

#####################################################################
# Nearest Neighbor Rank List                                        #
#####################################################################

@inbounds function Nearest_Neighbor_List{T<:Real}(Cost::Array{T,2},Pos::Int64,Num::Int64)
	new_cost = copy(Cost[Pos,:])
	tmp = maximum(new_cost)+1
	new_cost[Pos] = tmp
	NN_List::Vector{Int64} = Int64[]
	for i in 1:Num
		if minimum(new_cost) == tmp
			break
		end
		ind = findin(new_cost,minimum(new_cost))	
		push!(NN_List,ind...)
		new_cost[ind] = tmp
	end
	NN_List
end

@inbounds function Nearest_Neighbor_Rank_List{T<:Real}(Cost::Array{T,2})
	n::Int64 = size(Cost)[1]
	rank_list::Array{Int64,2} = zeros(Int64,n,n)
	for i in 1:n
		tmp = Nearest_Neighbor_List(Cost,i,n)
		for j in 1:length(tmp)
			rank_list[tmp[j],i] = j
		end
	end
	rank_list'
end

@inbounds function Nearest_Neighbor_Rank_List(data::BikeShareData)
	Nearest_Neighbor_Rank_List(data.Time)
end

#####################################################################
# Computing cost of a Tour                                          #
#####################################################################

@inbounds function time_of_tour(data::BikeShareData,Tour::Vector{Int64})
	Time::Float64 = 0.0
	for i in 1:length(Tour)-1
		Time += data.Time[Tour[i],Tour[i+1]]
	end
	Time
end

@inbounds function time_of_tour(data::BikeShareData,Tour::Vector{Int64},Instructions::Vector{Int64})
	Time::Float64 = float(sumabs(Instructions[2:end])*data.L_UL_time)
	for i in 1:length(Tour)-1
		Time += data.Time[Tour[i],Tour[i+1]]
	end
	Time
end

#####################################################################
# Various ranks in tour                                             #
#####################################################################

#####################################################################
# Edge rank list of a tour                                          #
#####################################################################

@inbounds function get_edge_rank_list(solution::BikeShareSolution)
	edge_rank_list::Vector{Vector{Int64}} = Vector{Int64}[]
	for v in 1:length(solution.Tour)
		push!(edge_rank_list,get_edge_rank_list(solution.Tour[v],solution.Rank_List))
	end
	edge_rank_list
end

@inbounds function get_edge_rank_list(tour::Vector{Int64},rank_list::Array{Int64,2})
	edge_rank_list::Vector{Int64} = Int64[]
	for i in 1:length(tour)-1
		push!(edge_rank_list,rank_list[tour[i],tour[i+1]])
	end
	edge_rank_list
end

#####################################################################
# Sorted edge rank list in a tour                                   #
#####################################################################

@inbounds function sorted_edge_rank_list(data::BikeShareData,solution::BikeShareSolution)
	serl::Vector{Vector{Int64}}, count::Vector{Int64} = Vector{Int64}[], Int64[]
	for v in 1:data.V
		tmp1,tmp2 = sorted_edge_rank_list(solution.Tour[v],solution.Rank_List,data.Time)
		push!(serl,tmp1); push!(count,tmp2)
	end
	serl,count
end

@inbounds function sorted_edge_rank_list(tour::Vector{Int64},rank_list::Array{Int64,2},cost_matrix::Array{Int64,2})
	edge_rank_list::Vector{Int64} = get_edge_rank_list(tour,rank_list)
	good_egde_rank::Float64 = mean(edge_rank_list) 
	count::Int64 = 0
	for i in edge_rank_list
		if i > good_egde_rank
			count += 1
		end
	end
	edge_cost_list::Vector{Int64}, edge_cost_rank_list::Vector{Int64}, serl::Vector{Int64} = Int64[], Int64[], Int64[]
	for i in 1:length(tour)-1
		push!(edge_cost_list,cost_matrix[tour[i],tour[i+1]])
	end
	for i in 1:length(edge_cost_list)
		push!(edge_cost_rank_list,indmax(edge_cost_list))
		edge_cost_list[indmax(edge_cost_list)] = 0
	end
	while countnz(edge_rank_list) > 0
		ind = findin(edge_rank_list,maximum(edge_rank_list))
		if length(ind) == 1
			push!(serl,ind[1])
			edge_rank_list[ind[1]] = 0
		else
			for ecrl in edge_cost_rank_list
				if ecrl in ind
					push!(serl,ecrl)
				end
			end
			edge_rank_list[ind] = 0
		end
	end
	serl,count
end

#####################################################################
# Highest edge rank in a tour                                       #
#####################################################################

@inbounds function highest_edge_rank(solution::BikeShareSolution)
	highest_rank::Vector{Int64}, edge_rank_list::Vector{Vector{Int64}} = Int64[], Vector{Int64}[]
	for v in 1:length(solution.Tour)
		push!(highest_rank,maximum(edge_rank_list[v]))
	end
	highest_rank
end

@inbounds function highest_edge_rank(Tour::Vector{Int64},Rank_List::Array{Int64,2})
	maximum(get_edge_rank_list(Tour,Rank_List))
end

#####################################################################
# Average edge rank in a tour                                       #
#####################################################################

@inbounds function average_edge_rank(solution::BikeShareSolution)
	average_rank::Vector{Float64}, edge_rank_list::Vector{Vector{Int64}} = Float64[], get_edge_rank_list(solution)
	for v in 1:length(solution.Tour)
		push!(average_rank,mean(edge_rank_list[v]))
	end
	average_rank
end

#####################################################################
# Node rank list of a tour                                          #
#####################################################################

@inbounds function get_node_rank_list(solution::BikeShareSolution)
	node_rank_list::Vector{Vector{Float64}} = Vector{Float64}[]
	for v in 1:length(solution.Tour)
		push!(node_rank_list,get_node_rank_list(solution.Tour[v],solution.Rank_List))
	end
	node_rank_list
end

@inbounds function get_node_rank_list(Tour::Vector{Int64},Rank_List::Array{Int64,2})
	node_rank_list::Vector{Float64} = zeros(length(Tour))
	for i in 2:length(Tour)-1
		node_rank_list[i] = (Rank_List[Tour[i-1],Tour[i]]+Rank_List[Tour[i],Tour[i+1]])/2
	end
	node_rank_list
end

#####################################################################
# Sorted edge rank list in a tour                                   #
#####################################################################

@inbounds function sorted_node_rank_list(data::BikeShareData,solution::BikeShareSolution)
	snrl::Vector{Vector{Int64}}, count::Vector{Int64} = Vector{Int64}[], Int64[]
	for v in 1:data.V
		tmp1,tmp2 = sorted_node_rank_list(solution.Tour[v],solution.Rank_List,data.Time)
		push!(snrl,tmp1); push!(count,tmp2)
	end
	snrl,count
end

@inbounds function sorted_node_rank_list(tour::Vector{Int64},rank_list::Array{Int64,2},cost_matrix::Array{Float64,2})
	node_rank_list::Vector{Float64} = get_node_rank_list(tour,rank_list)
	good_node_rank::Float64 = mean(node_rank_list[2:end-1])
	count::Int64 = 0
	for i in node_rank_list
		if i > good_node_rank
			count += 1
		end
	end
	node_cost_list::Vector{Float64}, node_cost_rank_list::Vector{Int64}, snrl::Vector{Int64} = Float64[], Int64[], Int64[]
	for i in 2:length(tour)-1
		push!(node_cost_list,cost_matrix[tour[i-1],tour[i]]+cost_matrix[tour[i],tour[i+1]])
	end
	for i in 1:length(node_cost_list)
		push!(node_cost_rank_list,indmax(node_cost_list))
		node_cost_list[indmax(node_cost_list)] = 0
	end
	node_cost_rank_list += 1
	while countnz(node_rank_list) > 0
		ind = findin(node_rank_list,maximum(node_rank_list))
		if length(ind) == 1
			push!(snrl,ind[1])
			node_rank_list[ind[1]] = 0
		else
			for ncrl in node_cost_rank_list
				if ncrl in ind
					push!(snrl,ncrl)
				end
			end
			node_rank_list[ind] = 0
		end
	end
	snrl,count
end

#####################################################################
# Highest node rank in a tour                                       #
#####################################################################

@inbounds function highest_node_rank(solution::BikeShareSolution)
	highest_rank::Vector{Int64} = Int64[]
	node_rank_list::Vector{Vector{Float64}} = get_node_rank_list(solution)
	for v in 1:length(solution.Tour)
		push!(highest_rank,maximum(node_rank_list[v]))
	end
	highest_rank
end

@inbounds function highest_node_rank(Tour::Vector{Int64},Rank_List::Array{Int64,2})
	maximum(get_node_rank_list(Tour,Rank_List))
end

#####################################################################
# Checking feasibility of a solution                                #
#####################################################################

@inbounds function check_feasibility(data::BikeShareData,solution::BikeShareSolution)
	O::Vector{Int64} = zeros(Int64,data.N)	
	for v in 1:data.V
		for i in 1:length(solution.Tour[v])-1
			if solution.Tour[v][i] == solution.Tour[v][i+1]
				println("Tour of vehicle $v is infeasible at $i")
				return "Current solution is Infeasible"
			end
		end
		q = 0
		for i in 1:length(solution.Instructions[v])
			q += solution.Instructions[v][i]
			if q > data.Q || q < 0
				println("Instruction of vehicle $v is infeasible at $i")
				return "Current solution is Infeasible"
    		end
    	end
    	for i in 2:length(solution.Instructions[v])-1
    	    O[solution.Tour[v][i]] += solution.Instructions[v][i]
    	end
    end
    for i in 1:data.N
    	if O[i] != data.O[i]
    		println("O is infeasible at $i")
        	return "Current solution is Infeasible"
        end
    end
    return "Current solution is Feasible"
end
