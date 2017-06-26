#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

@inbounds function optimal_bikes(cur::Int64,Q::Int64,O::Vector{Int64},node::Int64)
    n::Int64 = length(O)
    bikes::Vector{Int64} = fill(cur,n)
    for i in 1:n
    	bikes[i] += O[i]
    	if bikes[i] > Q
	    	bikes[i] = Q
		end
		if bikes[i] < 0
	    	bikes[i] = 0
		end
		bikes[i] -= cur
    end
    bikes[node] = 0
    bikes
end

#####################################################################
# Nearest Neighbor 1                                                #
#####################################################################

@inbounds function Nearest_Neighbor_1(node::Int64,Time::Array{Float64,2},OptB::Vector{Int64},L_t::Int64,O::Vector{Int64})
	alpha::Float64, ind::Vector{Int64} = 0.99, Int64[]
    if length(find(OptB)) == 0
		println("Optimal Bikes is a Zero Matrix")
		return 0
    end
    for i in 1:length(OptB)
    	if OptB[i] != 0
    		push!(ind,i)
   		end
   	end 
    ind[indmin(Time[node,ind])]
end

@inbounds function Nearest_Neighbor_1(node::Int64,Time::Array{Float64,2},OptB::Vector{Int64},L_t::Int64,O::Vector{Int64},rank_list::Array{Int64,2})
	if length(find(OptB)) == 0
		println("Optimal Bikes is a Zero Matrix")
		return 0
    end
    rank::Int64, ind::Int64 = size(Time)[1], 0
    for i in 1:length(OptB)
    	if OptB[i] != 0 && rank > rank_list[node,i]
    		ind = i
    		rank = rank_list[node,i]
   		end
   	end
   	ind
end

#####################################################################
# Nearest Neighbor 2                                                #
#####################################################################

@inbounds function Nearest_Neighbor_2(node::Int64,Time::Array{Float64,2},OptB::Vector{Int64},L_t::Int64,O::Vector{Int64})
    if length(find(OptB)) == 0
		println("Optimal Bikes is a Zero Matrix")
		return 0
    end
   	ind::Vector{Int64},cost::Vector{Float64} = Int64[], Float64[]
    for i in 1:length(OptB)
    	if OptB[i] != 0
    		push!(ind,i)
   		end
   	end
    for i in ind
    	push!(cost,abs(OptB[i])/Time[node,i])
    end
    ind[indmax(cost)]
end

@inbounds function Nearest_Neighbor_2(node::Int64,Time::Array{Float64,2},OptB::Vector{Int64},L_t::Int64,O::Vector{Int64},rank_list::Array{Int64,2})
    if length(find(OptB)) == 0
		println("Optimal Bikes is a Zero Matrix")
		return 0
    end
   ind::Vector{Int64},cost::Vector{Float64} = Int64[], Float64[]
    for i in 1:length(OptB)
    	if OptB[i] != 0
    		push!(ind,i)
   		end
   	end
    for i in ind
    	push!(cost,abs(OptB[i])/Time[node,i])
    end
    ind[indmax(cost)]
end

#####################################################################
# Nearest Neighbor 3                                                #
#####################################################################

@inbounds function Nearest_Neighbor_3(node::Int64,Time::Array{Float64,2},OptB::Vector{Int64},L_t::Int64,O::Vector{Int64})
	if length(find(OptB)) == 0
		println("Optimal Bikes is a Zero Matrix")
		return 0
    end
    ind::Vector{Int64} = Int64[]
    for i in 1:length(OptB)
    	if OptB[i] != 0
    		push!(ind,i)
   		end
   	end
    ind[rand(1:length(ind))]
end

@inbounds function Nearest_Neighbor_3(node::Int64,Time::Array{Float64,2},OptB::Vector{Int64},L_t::Int64,O::Vector{Int64},rank_list::Array{Int64,2})
	if length(find(OptB)) == 0
		println("Optimal Bikes is a Zero Matrix")
		return 0
    end
    ind::Vector{Int64} = Int64[]
    for i in 1:length(OptB)
    	if OptB[i] != 0
    		push!(ind,i)
   		end
   	end
    ind[rand(1:length(ind))]
end

#####################################################################
# Greedy tour construction heuristic                                #
#####################################################################

@inbounds function GCH(data::BikeShareData,Nearest_Neighbor::Symbol,starting_node::Int64,rank_list::Array{Int64,2})
	if starting_node == 0
		starting_node = rand(2:data.N)
    end
    
    tour::Vector{Int64}, instructions::Vector{Int64}, q::Vector{Int64}, O::Vector{Int64} = [1], [0], [0], copy(data.O)
    opt_bike::Vector{Int64}, nn::Int64 = optimal_bikes(q[end],data.Q,O,tour[end]), starting_node
    
    push!(tour,nn); push!(instructions,opt_bike[nn]); push!(q,q[end]+instructions[end]); O[nn] -= instructions[end]
			
    while length(find(O)) > 0
    	opt_bike = optimal_bikes(q[end],data.Q,O,tour[end])
		nn = eval(Nearest_Neighbor)(tour[end],data.Time,opt_bike,data.L_UL_time,O,rank_list)
		if nn == 0
	    	println("Somethings not right")
			break 	
		end
		push!(tour,nn); push!(instructions,opt_bike[nn]); push!(q,q[end]+instructions[end]); O[nn] -= instructions[end]
    end
    push!(tour,1); push!(instructions,0); push!(q,0)
	Time::Float64 = time_of_tour(data,tour,instructions)
	
	Tour::Vector{Vector{Int64}}, Instructions::Vector{Vector{Int64}} = Vector{Int64}[], Vector{Int64}[]
    push!(Tour,tour); push!(Instructions,instructions)
    
   	solution = BikeShareSolution(Time,[Time],Tour,Instructions,rank_list,0.0)
	ADJUST_INSTRUCTIONS(data,solution)
	solution	
end

#####################################################################
# Creating a Feasible Solution for Multiple vehicles                #
#####################################################################

#####################################################################
# Partition first rebalance second strategy                         #
#####################################################################

@inbounds function PFRBS(data::BikeShareData)
	t0::Float64 = time()
	
	initial_partitions::Vector{BikeShareData}, solutions::Vector{BikeShareSolution}, Rebalancing_Times::Vector{Float64}, Tour::Vector{Vector{Int64}}, Instructions::Vector{Vector{Int64}} = PBORMN(data), BikeShareSolution[], Float64[], Vector{Int64}[], Vector{Int64}[]
	for v in 1:data.V
		push!(solutions,INITIAL_SOLUTION_CREATION(initial_partitions[v]))
		push!(Rebalancing_Times,solutions[v].Makespan); push!(Tour,solutions[v].Tour[1]); push!(Instructions,solutions[v].Instructions[1])
	end
	Makespan::Float64 = maximum(Rebalancing_Times)
	
	t1::Float64 = time()
    
    BikeShareSolution(Makespan,Rebalancing_Times,Tour,Instructions,solutions[1].Rank_List,t1-t0)
end

#####################################################################
# Partitioning randomly                                             #
#####################################################################

@inbounds function PBORMN(data::BikeShareData)
	operations::Int64 = sum(abs(data.O))
	avg_operations::Int64, ind::Array{Int64,2}, Inv::Array{Int64,2}, Target::Array{Int64,2}, O::Array{Int64,2}, SI::Vector{Int64}, pos::Vector{Int64}, copy_O::Vector{Int64}, p::Int64, n::Int64, tmp::Int64 = fld(operations,data.V), zeros(Int64,data.N,data.V), zeros(Int64,data.N,data.V), zeros(Int64,data.N,data.V), zeros(Int64,data.N,data.V), zeros(Int64,data.V), shuffle(collect(2:data.N)), copy(data.O), 0, 0, 0
	
	if avg_operations%2 != 0
		avg_operations -= 1
	end
	 
	for v in 1:data.V-1
		p, n = round(Int64,avg_operations/2), -1*round(Int64,avg_operations/2)
		for i in 1:length(pos)
			if copy_O[pos[i]] > 0 && p > 0
				tmp = minimum([p,copy_O[pos[i]]])
				ind[pos[i],v], O[pos[i],v], Inv[pos[i],v] = 1, tmp, tmp
				p -= tmp
				copy_O[pos[i]] -= tmp
			else
				if copy_O[pos[i]] < 0 && n < 0
					tmp = maximum([n,copy_O[pos[i]]])
					ind[pos[i],v], O[pos[i],v], Target[pos[i],v] = 1, tmp, -1*tmp
					n -= tmp
					copy_O[pos[i]] -= tmp
				else
					if p == 0 && n == 0
						break
					end
				end
			end
		end
	end
	O[:,data.V] = copy_O[:]
	solution::Vector{BikeShareData} = BikeShareData[]
	for v in 1:data.V
		push!(solution,BikeShareData(data.N,O[:,v],1,data.Q,data.L_UL_time,data.Time))
	end
	solution
end

#####################################################################
# Repairing Infeasible Solutions                                    #
#####################################################################

#####################################################################
# Single Vehicle                                                    #
#####################################################################

@inbounds function TOUR_REPAIR(data::BikeShareData,Tour::Vector{Int64},Nearest_Neighbor::Symbol=:Nearest_Neighbor_1)
	O::Vector{Int64}, q::Int64, Instructions::Vector{Int64}, i::Int64 = copy(data.O), 0, [0], 2
	if Tour[1] != 1
		insert!(Tour,1,1)
	end
	while i <= length(Tour)
		opt_bike = optimal_bikes(q,data.Q,O,Tour[i-1])
		if opt_bike[Tour[i]] == 0
			deleteat!(Tour,i)
		else
			push!(Instructions,opt_bike[Tour[i]])
			q += Instructions[end]
			O[Tour[i]] -= opt_bike[Tour[i]]
			i += 1
		end
	end
	q = sum(Instructions)
	while sumabs(O) != 0
		opt_bike = optimal_bikes(q,data.Q,O,Tour[end])
		nn = eval(Nearest_Neighbor)(Tour[end],data.Time,opt_bike,data.L_UL_time,O)
		if nn == 0
	    	println("Somethings not right")
			break 	
		end
		push!(Tour,nn)
		push!(Instructions,opt_bike[nn])
		q += Instructions[end]
		O[nn] -= Instructions[end]
	end
	push!(Tour,1)
	push!(Instructions,0)
	(Tour,Instructions)
end

####################################################################
# Initial Solution Creation                                        #
####################################################################

@inbounds function INITIAL_SOLUTION_CREATION(data::BikeShareData)
	solution::BikeShareSolution = @match data.V begin
		1 => GCH(data,:Nearest_Neighbor_1,0,Nearest_Neighbor_Rank_List(data))
		_ => PFRBS(data)
	end
	VND(data,solution)
	solution
end
