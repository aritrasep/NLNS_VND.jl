#########################################################################################################
#  This file is part of the julia module for  computing rebalancing operations in bike sharing systems 	#
#  (c) Copyright 2015-17 by Aritra Pal                                                                  #
#  Permission is granted for academic research use.                                                     #
#  For other uses,  contact the author for licensing options.                                           #
#  Use at your own risk. I make no guarantees about the correctness or  usefulness of this code.        #
#########################################################################################################

module NLNS_VND

using Match
	
include("Types.jl")
include("Utilities.jl")
include("./Algorithm/SCRP/Overall_Algorithm.jl")
	
export BikeShareData, BikeShareSolution

export nlns_vnd

end
