#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

#####################################################################
# Delete-reinsert heuristic for tours of single vehicles            #
#####################################################################

@inbounds function INTRA_DELETE_REINSERT_LOCATIONS(data::BikeShareData,tour::Vector{Int64},instructions::Vector{Int64},q::Vector{Int64},pos::Int64)
	ind::Vector{Int64} = Int64[]
	for i in pos+1:length(tour)-1
		if instructions[pos] > 0 && q[pos+1] < instructions[pos]
			break
		end
		if instructions[pos] < 0 && q[pos+1] > data.Q + instructions[pos]
			break
		end
		if instructions[pos] > 0 && q[i] < instructions[pos]
			break
		end
		if instructions[pos] < 0 && q[i] > (data.Q + instructions[pos])
			break
		end
		push!(ind,i)
	end
	i = pos-1
	while i > 1
		if instructions[pos] > 0 && q[pos] > data.Q - instructions[pos]
			break
		end
		if instructions[pos] < 0 && q[pos] < -instructions[pos]
			break
		end
		if instructions[pos] > 0 && (q[i] > data.Q - instructions[pos] || q[i-1] > data.Q - instructions[pos])
			break
		end
		if instructions[pos] < 0 && (q[i] < -instructions[pos] || q[i-1] < -instructions[pos])
			break
		end
		push!(ind,i)
		i -= 1
	end
	sort!(ind)
end

@inbounds function INTRA_DELETE_REINSERT(data::BikeShareData,tour::Vector{Int64},instructions::Vector{Int64},rank_list::Array{Int64,2})

	q::Vector{Int64} = [instructions[1]]
	for i in 2:length(instructions)-1
		push!(q,instructions[i]+q[end])
	end	
	Time::Float64, highest_rank::Float64, pos::Int64, ind::Vector{Int64}, tmp::Float64 = time_of_tour(data,tour), highest_node_rank(tour,rank_list), 0, Int64[], 0.0	
	improvement_list::Vector{Float64} = [Time]
	
	while true
		i = 2
		while i < length(tour)
			if rank_list[tour[i-1],tour[i+1]] > highest_rank
				i += 1
				continue
			end
			Time, pos, ind = 0, i, INTRA_DELETE_REINSERT_LOCATIONS(data,tour,instructions,q,i)
			if length(ind) > 0
				for j in ind
					if j > i && minimum([rank_list[tour[j],tour[i]],rank_list[tour[i],tour[j+1]]]) > highest_rank
						continue
					end
					if j < i && minimum([rank_list[tour[j-1],tour[i]],rank_list[tour[i],tour[j]]]) > highest_rank
						continue
					end
					if j > i
						tmp = data.Time[tour[i-1],tour[i]] + data.Time[tour[i],tour[i+1]] + data.Time[tour[j],tour[j+1]] - data.Time[tour[i-1],tour[i+1]] - data.Time[tour[j],tour[i]] - data.Time[tour[i],tour[j+1]]
					end
					if j < i
						tmp = data.Time[tour[i-1],tour[i]] + data.Time[tour[i],tour[i+1]] + data.Time[tour[j],tour[j-1]] - data.Time[tour[i-1],tour[i+1]] - data.Time[tour[j-1],tour[i]] - data.Time[tour[i],tour[j]]
					end
					if tmp > Time
						Time, pos = tmp, j
					end
				end	
				if pos != i
					if pos > i
						tour[i:pos-1], tour[pos], instructions[i:pos-1], instructions[pos] = tour[i+1:pos], tour[i], instructions[i+1:pos], instructions[i]
					else
						tour[pos+1:i], tour[pos], instructions[pos+1:i], instructions[pos] = tour[pos:i-1], tour[i], instructions[pos:i-1], instructions[i]	
					end
					k = 2
					while k < length(tour)
						if tour[k] == tour[k+1]
							instructions[k] += instructions[k+1]
							deleteat!(tour,k+1); deleteat!(instructions,k+1)
						end
						k += 1
					end
					q = [instructions[1]]
					for k in 2:length(instructions)-1
						push!(q,instructions[k]+q[end])
					end
					highest_rank = highest_node_rank(tour,rank_list)
				else
					i += 1
				end
			else
				i += 1
			end
		end
		Time = time_of_tour(data,tour)
		push!(improvement_list,Time)
		if length(improvement_list) > 2
			if improvement_list[end] == improvement_list[end-1]
				break
			end
		end
	end
	(tour,instructions)
end

@inbounds function INTRA_DELETE_REINSERT(data::BikeShareData,solution::BikeShareSolution)
	time_::Float64 = 0
	for v in 1:data.V
		time_ = time_of_tour(data,solution.Tour[v],solution.Instructions[v])
		(solution.Tour[v],solution.Instructions[v]) = INTRA_DELETE_REINSERT(data,solution.Tour[v],solution.Instructions[v],solution.Rank_List)
		solution.Rebalancing_Times[v] = time_
	end
	solution.Makespan = maximum(solution.Rebalancing_Times)
end

#####################################################################
# Interchanging 2 sequences of edges                                #
#####################################################################

@inbounds function INTRA_ARC_EXCHANGE(data::BikeShareData,tour::Vector{Int64},instructions::Vector{Int64},rank_list::Array{Int64,2})
	
	q::Vector{Int64} = [instructions[1]]
	for i in 2:length(instructions)-1
		push!(q,instructions[i]+q[end])
	end	
	Time::Float64, highest_rank::Int64, pos1::Vector{Int64}, pos2::Vector{Int64}, Pos1::Int64, Pos2::Int64, tmp::Float64, tmp1::Float64, tmp2::Vector{Int64}, step::Int64  = time_of_tour(data,tour), highest_edge_rank(tour,rank_list), Int64[], Int64[], 0, 0, 0.0, 0.0, Int64[], 1
	
	improvement_list::Vector{Float64} = [Time]

	while (data.N >= 200 && step <= data.N/20) || (data.N < 200 && step <= length(q)-4)
		i = 1
		while i <= length(q)-2-step
			pos1, pos2, Pos1, Pos2, tmp = findin(q[i+step+1:end-1],q[i])+i+step, findin(q[i+step+2:end],q[i+step])+i+step+1, 0, 0, 0
			for p1 in pos1, p2 in pos2
				if p1 < p2
					if minimum([rank_list[tour[i],tour[p1+1]],rank_list[tour[p2],tour[i+step+1]],rank_list[tour[p1],tour[i+1]],rank_list[tour[i+step],tour[p2+1]]]) >= highest_rank
						continue
					end
					tmp1 = ( data.Time[tour[i],tour[i+1]] + data.Time[tour[i+step],tour[i+step+1]] + data.Time[tour[p1],tour[p1+1]] + data.Time[tour[p2],tour[p2+1]] ) - ( data.Time[tour[i],tour[p1+1]] + data.Time[tour[p2],tour[i+step+1]] + data.Time[tour[p1],tour[i+1]] + data.Time[tour[i+step],tour[p2+1]] )			
					if tmp1 > tmp
						tmp, Pos1, Pos2 = tmp1, p1, p2
					end
				end
			end
			if Pos1 > 0 && Pos2 > 0
				tmp2 = splice!(tour,Pos1+1:Pos2,tour[i+1:i+step]); splice!(tour,i+1:i+step,tmp2)
				tmp2 = splice!(instructions,Pos1+1:Pos2,instructions[i+1:i+step]); splice!(instructions,i+1:i+step,tmp2)
				k = 2
				while k < length(tour)
					if tour[k] == tour[k+1]
						instructions[k] += instructions[k+1]
						deleteat!(tour,k+1); deleteat!(instructions,k+1)
					end
					k += 1
				end
				q = [instructions[1]]
				for k in 2:length(instructions)-1
					push!(q,instructions[k]+q[end])
				end
				highest_rank = highest_edge_rank(tour,rank_list)
			end
			i += 1
		end
		Time = time_of_tour(data,tour)
		push!(improvement_list,Time)
		if improvement_list[end] == improvement_list[end-1]
			step += 1
		end
	end
	(tour,instructions)
end

@inbounds function INTRA_ARC_EXCHANGE(data::BikeShareData,solution::BikeShareSolution)
	time_::Float64 = 0
	for v in 1:data.V
		time_ = time_of_tour(data,solution.Tour[v],solution.Instructions[v])
		(solution.Tour[v],solution.Instructions[v]) = INTRA_ARC_EXCHANGE(data,solution.Tour[v],solution.Instructions[v],solution.Rank_List)
		solution.Rebalancing_Times[v] = time_
	end
	solution.Makespan = maximum(solution.Rebalancing_Times)
end

#####################################################################
# Intra-swapping / swapping nodes inside a tour of a vehicle        #
#####################################################################

@inbounds function INTRA_NODE_SWAP(data::BikeShareData,tour::Vector{Int64},instructions::Vector{Int64},rank_list::Array{Int64,2})
	
	q::Vector{Int64} = [instructions[1]]
	for i in 2:length(instructions)-1
		push!(q,instructions[i]+q[end])
	end	
	Time::Float64, highest_rank::Float64, pos::Vector{Int64}, pos2::Int64, tmp::Float64 = time_of_tour(data,tour), highest_node_rank(tour,rank_list), Int64[], 0, 0.0
	improvement_list::Vector{Float64} = [Time]
	while true
		for i in 2:length(tour)-2
			if i > length(instructions)
				break
			end
			j = instructions[i]::Int64
			pos, pos2, tmp = findin(instructions[i+1:end-1],j)+i, 0, 0
			for k in pos
				if tour[k] != tour[i]
					if minimum([rank_list[tour[i-1],tour[k]],rank_list[tour[k],tour[i+1]],rank_list[tour[k-1],tour[i]],rank_list[tour[i],tour[k+1]]]) >= highest_rank
						continue
					end
					tour[i],tour[k] = tour[k],tour[i]
					tmp1 = time_of_tour(data,tour,instructions)
					tour[i],tour[k] = tour[k],tour[i]
					if Time-tmp1 > tmp
						tmp, pos2 = Time-tmp1, k
					end
				end
			end
			if tmp > 0
				tour[i],tour[pos2] = tour[pos2],tour[i]
				k = 2
				while k < length(tour)
					if tour[k] == tour[k+1]
						instructions[k] += instructions[k+1]
						deleteat!(tour,k+1); deleteat!(instructions,k+1)
					end
					k += 1
				end
				q = [instructions[1]]
				for k in 2:length(instructions)-1
					push!(q,instructions[k]+q[end])
				end
				highest_rank = highest_node_rank(tour,rank_list)
			end
		end
		Time = time_of_tour(data,tour)
		push!(improvement_list,Time)
		if length(improvement_list) > 2
			if improvement_list[end] == improvement_list[end-1] 
				break
			end
		end
	end
	(tour,instructions)
end

@inbounds function INTRA_NODE_SWAP(data::BikeShareData,solution::BikeShareSolution)
	time_::Float64 = 0
	for v in 1:data.V
		time_ = time_of_tour(data,solution.Tour[v],solution.Instructions[v])
		(solution.Tour[v],solution.Instructions[v]) = INTRA_NODE_SWAP(data,solution.Tour[v],solution.Instructions[v],solution.Rank_List)
		solution.Rebalancing_Times[v] = time_
	end
	solution.Makespan = maximum(solution.Rebalancing_Times)
end

@inbounds function INTRA_NODE_SWAP2(data::BikeShareData,tour::Vector{Int64},instructions::Vector{Int64},rank_list::Array{Int64,2})
	
	q::Vector{Int64} = [instructions[1]]
	for i in 2:length(instructions)-1
		push!(q,instructions[i]+q[end])
	end	
	Time::Float64, highest_rank::Float64, pos::Int64, tmp::Float64, tmp_time::Float64 = time_of_tour(data,tour), highest_node_rank(tour,rank_list), 0, 0.0, 0.0
	improvement_list::Vector{Float64} = [Time]
	while true
		for i in 2:length(tour)-2
			pos, tmp = 0, 0
			for j in i+2:length(tour)-1
				if tour[i] != tour[j] && tour[i] != tour[j-1] && tour[i] != tour[j+1] && tour[j] != tour[i-1] && tour[j] != tour[i+1]
					if minimum([rank_list[tour[i-1],tour[j]],rank_list[tour[j],tour[i+1]],rank_list[tour[j-1],tour[i]],rank_list[tour[i],tour[j+1]]]) >= highest_rank
						continue
					end
					tmp_time = 0
					if instructions[j]-instructions[i] > 0 && data.Q-maximum(q[i:j-1]) >= instructions[j]-instructions[i]
						tmp_time = ( data.Time[tour[i-1],tour[i]] + data.Time[tour[i],tour[i+1]] + data.Time[tour[j-1],tour[j]] + data.Time[tour[j],tour[j+1]] ) - ( data.Time[tour[i-1],tour[j]] + data.Time[tour[j],tour[i+1]] + data.Time[tour[j-1],tour[i]] + data.Time[tour[i],tour[j+1]] )
					elseif instructions[j]-instructions[i] < 0 && minimum(q[i:j-1]) >= instructions[i]-instructions[j]
						tmp_time = ( data.Time[tour[i-1],tour[i]] + data.Time[tour[i],tour[i+1]] + data.Time[tour[j-1],tour[j]] + data.Time[tour[j],tour[j+1]] ) - ( data.Time[tour[i-1],tour[j]] + data.Time[tour[j],tour[i+1]] + data.Time[tour[j-1],tour[i]] + data.Time[tour[i],tour[j+1]] )
					end
					if tmp_time > tmp
						tmp, pos = tmp_time, j
					end
				end
			end
			if tmp > 0
				tour[i], tour[pos], instructions[i], instructions[pos] = tour[pos], tour[i], instructions[pos], instructions[i]
				k = 2
				while k < length(tour)
					if tour[k] == tour[k+1]
						instructions[k] += instructions[k+1]
						deleteat!(tour,k+1); deleteat!(instructions,k+1)
					end
					k += 1
				end
				q = [instructions[1]]
				for k in 2:length(instructions)-1
					push!(q,instructions[k]+q[end])
				end
				highest_rank = highest_node_rank(tour,rank_list)
			end
		end
		Time = time_of_tour(data,tour)
		push!(improvement_list,Time)
		if length(improvement_list) > 2
			if improvement_list[end] == improvement_list[end-1] 
				break
			end
		end
	end
	(tour,instructions)
end

@inbounds function INTRA_NODE_SWAP2(data::BikeShareData,solution::BikeShareSolution)
	time_::Float64 = 0
	for v in 1:data.V
		time_ = time_of_tour(data,solution.Tour[v],solution.Instructions[v])
		(solution.Tour[v],solution.Instructions[v]) = INTRA_NODE_SWAP2(data,solution.Tour[v],solution.Instructions[v],solution.Rank_List)
		solution.Rebalancing_Times[v] = time_
	end
	solution.Makespan = maximum(solution.Rebalancing_Times)
end

#####################################################################
# Adjusting instructions of nodes in the tour                       #
#####################################################################

@inbounds function ADJUST_INSTRUCTIONS(data::BikeShareData,tour::Vector{Int64},instructions::Vector{Int64})
	
	q::Vector{Int64} = [instructions[1]]
	for i in 2:length(instructions)-1
		push!(q,instructions[i]+q[end])
	end
	pos::Dict{Int64,Vector{Int64}}, tmp::Vector{Int64}, tmp1::Int64, tmpI1::Int64, tmp2::Int64, tmpI2::Int64, UP::Int64, DOWN::Int64 = Dict{Int64,Vector{Int64}}(), Int64[], 0, 0, 0, 0, 0, 0
	for i in 2:length(tour)-3
		tmp = findin(tour,tour[i])
		if length(tmp) > 1
			pos[tour[i]] = tmp
		end
	end
	for i in collect(keys(pos))
		for j in 1:length(pos[i])-1
			tmp1 = pos[i][j]
			tmpI1 = instructions[tmp1]  
			for k in j+1:length(pos[i])
				tmp2 = pos[i][k]
				tmpI2, min_q, max_q = instructions[tmp2], minimum(q[tmp1:tmp2-1]), maximum(q[tmp1:tmp2-1])
				if tmpI1*tmpI2 >= 0
					UP, DOWN = 0, 0
					if tmpI1 > 0
						UP, DOWN = minimum([data.Q-tmpI1,tmpI2,data.Q-max_q]), minimum([tmpI1,data.Q-tmpI2,min_q])
						if UP > DOWN
							instructions[tmp1] += UP
							instructions[tmp2] -= UP
							q = [instructions[1]]
							for l in 2:length(instructions)-1
								push!(q,instructions[l]+q[end])
							end
						else
							if DOWN > UP
								instructions[tmp1] -= DOWN
								instructions[tmp2] += DOWN
								q = [instructions[1]]
								for l in 2:length(instructions)-1
									push!(q,instructions[l]+q[end])
								end
							end
						end
					end
					if tmpI1 < 0
						UP, DOWN = -1*minimum([abs(tmpI2),min_q]), -1*minimum([abs(tmpI1),data.Q-max_q])
						if UP < DOWN
							instructions[tmp1] += UP
							instructions[tmp2] -= UP
							q = [instructions[1]]
							for l in 2:length(instructions)-1
								push!(q,instructions[l]+q[end])
							end
						else
							if DOWN < UP
								instructions[tmp1] -= DOWN
								instructions[tmp2] += DOWN
								q = [instructions[1]]
								for l in 2:length(instructions)-1
									push!(q,instructions[l]+q[end])
								end
							end
						end	
					end
				end
			end
		end
	end
	tmp_pos::Vector{Int64}, ind::Vector{Int64} = findin(instructions[2:end-1],0)+1, Int64[]
	if length(tmp_pos) > 0
		for tp in tmp_pos
			if data.Time[tour[tp-1],tour[tp]]+data.Time[tour[tp],tour[tp+1]] <= data.Time[tour[tp-1],tour[tp+1]]
				push!(ind,tp)
			end
		end
		if length(ind) > 0
			tmp_pos = sort(setdiff(tmp_pos,ind))
		end
	end
	deleteat!(tour,tmp_pos); deleteat!(instructions,tmp_pos)
end	

@inbounds function ADJUST_INSTRUCTIONS(data::BikeShareData,solution::BikeShareSolution)
	time_::Float64 = 0
	for v in 1:data.V
		time_ = time_of_tour(data,solution.Tour[v],solution.Instructions[v])
		ADJUST_INSTRUCTIONS(data,solution.Tour[v],solution.Instructions[v])
		solution.Rebalancing_Times[v] = time_
	end
	solution.Makespan = maximum(solution.Rebalancing_Times)
end
