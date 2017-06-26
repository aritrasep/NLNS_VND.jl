#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

import Base.copy

#####################################################################
# Julia type for storing data for a test case                       #
#####################################################################

type BikeShareData
    N::Int64
    O::Vector{Int64}
	V::Int64
    Q::Int64
    L_UL_time::Int64
    Time::Array{Float64,2}
end

@inbounds function copy(instance::BikeShareData)
    BikeShareData(instance.N,instance.O,instance.V,instance.Q,instance.L_UL_time,instance.Time)
end

#####################################################################
# Julia type for storing solutions of a test case                   #
#####################################################################

type BikeShareSolution
    Makespan::Float64
    Rebalancing_Times::Vector{Float64}
    Tour::Vector{Vector{Int64}}
    Instructions::Vector{Vector{Int64}}
    Rank_List::Array{Int64,2}
    Computing_Time::Float64
end

@inbounds function copy(instance::BikeShareSolution)
    BikeShareSolution(instance.Makespan,instance.Rebalancing_Times,instance.Tour,instance.Instructions,instance.Rank_List,instance.condition,instance.Computing_Time)
end
