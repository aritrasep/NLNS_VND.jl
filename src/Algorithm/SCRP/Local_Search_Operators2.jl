#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

#####################################################################
# 2-OPT Neighborhood                                                #
#####################################################################

@inbounds function INTER_2_OPT(data::BikeShareData,Tour1::Vector{Int64},Instructions1::Vector{Int64},Tour2::Vector{Int64},Instructions2::Vector{Int64})
	
	q1::Vector{Int64} = [Instructions1[1]]
	for i in 2:length(Instructions1)-1
		push!(q1,Instructions1[i]+q1[end])
	end
	q2::Vector{Int64} = [Instructions2[1]]
	for i in 2:length(Instructions2)-1
		push!(q2,Instructions2[i]+q2[end])
	end
	Time1::Float64, Time2::Float64, pos::Int64, tmp::Vector{Int64}, tmp_time1::Float64, tmp_time2::Float64, tmp_makespan::Float64, improvement_list::Vector{Float64} = time_of_tour(data,Tour1,Instructions1), time_of_tour(data,Tour2,Instructions2), 0, Int64[], 0.0, 0.0, 0.0, Float64[]
	Makespan::Float64 = Time1>=Time2?Time1:Time2 				
	while true
		i = 1
		while i < length(Tour1)
			pos = 0
			locations = findin(q2,q1[i])
			for j in locations
				#####################################################
				# Condition for Granularity needs to be added here  #
				#####################################################
				tmp_tour1, tmp_instructions1, tmp_tour2, tmp_instructions2 = copy(Tour1), copy(Instructions1), copy(Tour2), copy(Instructions2)
				tmp = splice!(tmp_tour1,i+1:length(tmp_tour1),tmp_tour2[j+1:end])
				splice!(tmp_tour2,j+1:length(tmp_tour2),tmp)
				tmp = splice!(tmp_instructions1,i+1:length(tmp_instructions1),tmp_instructions2[j+1:end])
				splice!(tmp_instructions2,j+1:length(tmp_instructions2),tmp)
				tmp_time1, tmp_time2 = time_of_tour(data,tmp_tour1,tmp_instructions1), time_of_tour(data,tmp_tour2,tmp_instructions2)
				tmp_makespan = tmp_time1 >= tmp_time2?tmp_time1:tmp_time2
				if tmp_makespan < Makespan
					Makespan, pos = tmp_makespan, j
				end
			end
			if pos > 0
				tmp = splice!(Tour1,i+1:length(Tour1),Tour2[pos+1:end])
				splice!(Tour2,pos+1:length(Tour2),tmp)
				tmp = splice!(Instructions1,i+1:length(Instructions1),Instructions2[pos+1:end])
				splice!(Instructions2,pos+1:length(Instructions2),tmp)
				k = 2
				while k < length(Tour1)
					if Tour1[k] == Tour1[k+1]
						Instructions1[k] += Instructions1[k+1]
						deleteat!(Tour1,k+1); deleteat!(Instructions1,k+1)
					end
					k += 1
				end
				k = 2
				while k < length(Tour2)
					if Tour2[k] == Tour2[k+1]
						Instructions2[k] += Instructions2[k+1]
						deleteat!(Tour2,k+1); deleteat!(Instructions2,k+1)
					end
					k += 1
				end
				q1 = [Instructions1[1]]
				for i in 2:length(Instructions1)-1
					push!(q1,Instructions1[i]+q1[end])
				end
				q2 = [Instructions2[1]]
				for i in 2:length(Instructions2)-1
					push!(q2,Instructions2[i]+q2[end])
				end
				Time1 = time_of_tour(data,Tour1,Instructions1); Time2 = time_of_tour(data,Tour2,Instructions2)
				Makespan = Time1 >= Time2?Time1:Time2
			end
			i += 1
		end
		push!(improvement_list,Time1>=Time2?Time1:Time2)
		if length(improvement_list) > 1 && improvement_list[end] == improvement_list[end-1]
			break
		end
	end
	(Tour1,Instructions1,Time1,Tour2,Instructions2,Time2)
end

@inbounds function INTER_2_OPT(data::BikeShareData,solution::BikeShareSolution)
	improvement_list = [solution.Makespan]
	
	tmp_time1::Float64, tmp_time2::Float64 = 0.0, 0.0
	#println("Initial Makespan is $(solution.Makespan) INTER_2_OPT")
	while true
	 	for i in 1:data.V-1
	 		for j in i+1:data.V
	 			(solution.Tour[i],solution.Instructions[i],tmp_time1,solution.Tour[j],solution.Instructions[j],tmp_time2) = INTER_2_OPT(data,solution.Tour[i],solution.Instructions[i],solution.Tour[j],solution.Instructions[j])
	 			if tmp_time1 != solution.Rebalancing_Times[i] || tmp_time2 != solution.Rebalancing_Times[j]
	 				solution.Rebalancing_Times[i], solution.Rebalancing_Times[j] = tmp_time1, tmp_time2
	 				solution.Makespan = maximum(solution.Rebalancing_Times)
	 			end
	 		end
	 	end
	 	push!(improvement_list,solution.Makespan)
	 	if improvement_list[end] == improvement_list[end-1]
	 		break
	 	end
	 end
	 #println("Final Makespan is $(solution.Makespan) INTER_2_OPT")
	 #if check_feasibility(data,solution) != "Current solution is Feasible"
	 #	println("INTER_2_OPT fucked up")
	 #end
end	

#####################################################################
# Inter-swapping / swapping nodes amongst tours of vehicles         #
#####################################################################

@inbounds function INTER_NODE_SWAP(data::BikeShareData,Tour1::Vector{Int64},Instructions1::Vector{Int64},Tour2::Vector{Int64},Instructions2::Vector{Int64})
	Time1::Float64, Time2::Float64, pos::Int64, tmp_time1::Float64, tmp_time2::Float64, tmp_makespan::Float64, improvement_list::Vector{Float64} = time_of_tour(data,Tour1,Instructions1), time_of_tour(data,Tour2,Instructions2), 0, 0.0, 0.0, 0.0, Float64[]
	Makespan::Float64 = Time1>=Time2?Time1:Time2

	while true
		i = 2
		while i < length(Tour1)
			pos = 0
			j = 2
			while j < length(Tour2)
				if i >= length(Instructions1) || j >= length(Instructions2)
					break
				end
				if Instructions1[i] == Instructions2[j] && Tour1[i] != Tour2[j]
					#####################################################
					# Condition for Granularity needs to be added here  #
					#####################################################
					if data.Time[Tour1[i-1],Tour1[i]] + data.Time[Tour1[i],Tour1[i+1]] + data.Time[Tour2[j-1],Tour2[j]] + data.Time[Tour2[j],Tour2[j+1]] > data.Time[Tour1[i-1],Tour2[j]] + data.Time[Tour2[j],Tour1[i+1]] + data.Time[Tour2[j-1],Tour1[i]] + data.Time[Tour1[i],Tour2[j+1]]
						tmp_time1 = Time1 - data.Time[Tour1[i-1],Tour1[i]] - data.Time[Tour1[i],Tour1[i+1]] + data.Time[Tour1[i-1],Tour2[j]] + data.Time[Tour2[j],Tour1[i+1]]
						tmp_time2 = Time2 - data.Time[Tour2[j-1],Tour2[j]] - data.Time[Tour2[j],Tour2[j+1]] + data.Time[Tour2[j-1],Tour1[i]] + data.Time[Tour1[i],Tour2[j+1]]
						tmp_makespan = tmp_time1>=tmp_time2?tmp_time1:tmp_time2
						if Makespan > tmp_makespan
							Makespan, pos = tmp_makespan, j
						end
					end
				end
				j += 1
			end
			if pos > 0
				Tour1[i], Tour2[pos] = Tour2[pos], Tour1[i]
				if Tour1[i-1] == Tour1[i] || Tour1[i] == Tour1[i+1]
					k = 2
					while k < length(Tour1)
						if Tour1[k] == Tour1[k+1]
							Instructions1[k] += Instructions1[k+1]
							deleteat!(Tour1,k+1); deleteat!(Instructions1,k+1)
						end
						k += 1
					end
				end
				if Tour2[pos-1] == Tour2[pos] || Tour2[pos] == Tour2[pos+1]
					k = 2
					while k < length(Tour2)
						if Tour2[k] == Tour2[k+1]
							Instructions2[k] += Instructions2[k+1]
							deleteat!(Tour2,k+1); deleteat!(Instructions2,k+1)
						end
						k += 1
					end
				end
				Time1, Time2 = time_of_tour(data,Tour1,Instructions1), time_of_tour(data,Tour2,Instructions2)
			end
			i += 1
		end
		push!(improvement_list,Time1 >= Time2?Time1:Time2)
		if length(improvement_list) > 1 && improvement_list[end] == improvement_list[end-1]
			break
		end
	end
	(Tour1,Instructions1,Time1,Tour2,Instructions2,Time2)
end

@inbounds function INTER_NODE_SWAP(data::BikeShareData,solution::BikeShareSolution)
	improvement_list::Vector{Float64} = [solution.Makespan]
	
	tmp_time1::Float64, tmp_time2::Float64 = 0.0, 0.0
	#println("Initial Makespan is $(solution.Makespan) INTER_NODE_SWAP")
	while true
	 	for i in 1:data.V-1
	 		for j in i+1:data.V
	 			(solution.Tour[i],solution.Instructions[i],tmp_time1,solution.Tour[j],solution.Instructions[j],tmp_time2) = INTER_NODE_SWAP(data,solution.Tour[i],solution.Instructions[i],solution.Tour[j],solution.Instructions[j])
	 			if tmp_time1 != solution.Rebalancing_Times[i] || tmp_time2 != solution.Rebalancing_Times[j]
	 				solution.Rebalancing_Times[i], solution.Rebalancing_Times[j] = tmp_time1, tmp_time2
	 				solution.Makespan = maximum(solution.Rebalancing_Times)
	 			end
	 		end
	 	end
	 	push!(improvement_list,solution.Makespan)
	 	if improvement_list[end] == improvement_list[end-1]
	 		break
	 	end
	 end
	 #println("Final Makespan is $(solution.Makespan) INTER_NODE_SWAP")
	 #if check_feasibility(data,solution) != "Current solution is Feasible"
	 #	println("INTER_NODE_SWAP fucked up")
	 #end
end	

#####################################################################
# Exchanging sequences                                              #
#####################################################################

@inbounds function INTER_ARC_EXCHANGE(data::BikeShareData,Tour1::Vector{Int64},Instructions1::Vector{Int64},Tour2::Vector{Int64},Instructions2::Vector{Int64})
	q_initial::Int64, q_final::Int64, q_2::Vector{Int64}, Time1::Float64, Time2::Float64, pos::Int64, tmp::Vector{Int64}, tmp_time1::Float64, tmp_time2::Float64, tmp_makespan::Float64, improvement_list::Vector{Float64}, pos1::Vector{Int64}, pos2::Vector{Int64}, tmp_pos1::Int64, tmp_pos2::Int64, step_limit::Int64, step::Int64 = 0, 0, [Instructions2[1]], time_of_tour(data,Tour1,Instructions1), time_of_tour(data,Tour2,Instructions2), 0, Int64[], 0.0, 0.0, 0.0, Float64[], Int64[], Int64[], 0, 0, length(Tour1)-2, 1
	Makespan::Float64 = Time1>=Time2?Time1:Time2
	for i in 2:length(Instructions2)-1
		push!(q_2,q_2[end]+Instructions2[i])
	end
	while true
		i = 1
		while i <= length(Tour1)-1-step
			tmp_makespan = Time1>=Time2?Time1:Time2
			q_initial, q_final = sum(Instructions1[1:i]), sum(Instructions1[1:i+step])
			pos1 = findin(q_2,q_initial)
			if length(pos1) > 0
				tmp_pos1, tmp_pos2 = 0, 0
				for j in pos1						
					pos2 = findin(q_2[j+1:end],q_final)+j
					if length(pos2) > 0
						for k in pos2
							if k >= length(Tour2)
								break
							end
							#####################################################
							# Condition for Granularity needs to be added here  #
							#####################################################
							tmp = splice!(Tour1,i+1:i+step,Tour2[j+1:k]); splice!(Tour2,j+1:k,tmp)
							tmp = splice!(Instructions1,i+1:i+step,Instructions2[j+1:k]); splice!(Instructions2,j+1:k,tmp)
							tmp_time1, tmp_time2 = time_of_tour(data,Tour1,Instructions1), time_of_tour(data,Tour2,Instructions2)
							if tmp_makespan > (tmp_time1>=tmp_time2?tmp_time1:tmp_time2)
								tmp_makespan = tmp_time1>=tmp_time2?tmp_time1:tmp_time2
								tmp_pos1, tmp_pos2 = j, k
							end
							tmp = splice!(Tour1,i+1:i+k-j,Tour2[j+1:j+step]); splice!(Tour2,j+1:j+step,tmp)
							tmp = splice!(Instructions1,i+1:i+k-j,Instructions2[j+1:j+step]); splice!(Instructions2,j+1:j+step,tmp)
						end
					end
				end
				if tmp_pos1 > 0 && tmp_pos2 > 0
					tmp = splice!(Tour1,i+1:i+step,Tour2[tmp_pos1+1:tmp_pos2]); splice!(Tour2,tmp_pos1+1:tmp_pos2,tmp)
					tmp = splice!(Instructions1,i+1:i+step,Instructions2[tmp_pos1+1:tmp_pos2]); splice!(Instructions2,tmp_pos1+1:tmp_pos2,tmp)
					k = 2
					while k < length(Tour1)
						if Tour1[k] == Tour1[k+1]
							Instructions1[k] += Instructions1[k+1]
							deleteat!(Tour1,k+1); deleteat!(Instructions1,k+1)
						end
						k += 1
					end
					k = 2
					while k < length(Tour2)
						if Tour2[k] == Tour2[k+1]
							Instructions2[k] += Instructions2[k+1]
							deleteat!(Tour2,k+1); deleteat!(Instructions2,k+1)
						end
						k += 1
					end
					q_2 = [Instructions2[1]]
					for i in 2:length(Instructions2)-1
						push!(q_2,q_2[end]+Instructions2[i])
					end	
					Time1, Time2 = time_of_tour(data,Tour1,Instructions1), time_of_tour(data,Tour2,Instructions2)
				end
			end
			i += 1
		end
		step += 1
		push!(improvement_list,Time1 >= Time2?Time1:Time2)
		if length(improvement_list) > 1 && improvement_list[end] == improvement_list[end-1]
			break
		end
	end
	(Tour1,Instructions1,Time1,Tour2,Instructions2,Time2)
end

@inbounds function INTER_ARC_EXCHANGE(data::BikeShareData,solution::BikeShareSolution)
	improvement_list::Vector{Float64} = [solution.Makespan]
	
	tmp_time1::Float64, tmp_time2::Float64 = 0.0, 0.0
	#println("Initial Makespan is $(solution.Makespan) INTER_ARC_EXCHANGE")
	while true
	 	for i in 1:data.V-1
	 		for j in i+1:data.V
	 			(solution.Tour[i],solution.Instructions[i],tmp_time1,solution.Tour[j],solution.Instructions[j],tmp_time2) = INTER_ARC_EXCHANGE(data,solution.Tour[i],solution.Instructions[i],solution.Tour[j],solution.Instructions[j])
	 			if tmp_time1 != solution.Rebalancing_Times[i] || tmp_time2 != solution.Rebalancing_Times[j]
	 				solution.Rebalancing_Times[i], solution.Rebalancing_Times[j] = tmp_time1, tmp_time2
	 				solution.Makespan = maximum(solution.Rebalancing_Times)
	 			end
	 		end
	 	end
	 	push!(improvement_list,solution.Makespan)
	 	if improvement_list[end] == improvement_list[end-1]
	 		break
	 	end
	 end
	 #println("Final Makespan is $(solution.Makespan) INTER_ARC_EXCHANGE")
	 #if check_feasibility(data,solution) != "Current solution is Feasible"
	 #	println("INTER_ARC_EXCHANGE fucked up")	
	 #end
end
