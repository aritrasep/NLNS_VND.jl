#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

#####################################################################
# Complex Neighborhood Search / Variable Neighborhood Descent       #
#####################################################################

@inbounds function VND(data::BikeShareData,solution::BikeShareSolution)

	neighborhoods::Vector{Symbol} = @match data.V begin
		1 => [:INTRA_DELETE_REINSERT,:INTRA_NODE_SWAP,:INTRA_NODE_SWAP2,:INTRA_ARC_EXCHANGE,:ADJUST_INSTRUCTIONS]
		_ => [:INTRA_DELETE_REINSERT,:INTRA_NODE_SWAP,:INTRA_NODE_SWAP2,:INTRA_ARC_EXCHANGE,:INTER_2_OPT,:INTER_NODE_SWAP,:INTER_ARC_EXCHANGE,:ADJUST_INSTRUCTIONS]
	end
	VND(data,solution,neighborhoods)
	
end

@inbounds function VND(data::BikeShareData,solution::BikeShareSolution,neighborhoods::Vector{Symbol})
	improvement_list::Vector{Float64} = [solution.Makespan]
	while true
		for i in 1:length(neighborhoods)
			eval(neighborhoods[i])(data,solution)
		end
		push!(improvement_list,solution.Makespan)
		if improvement_list[end] == improvement_list[end-1]
			break
		end
	end
end
