export Population, PopList, ConfigInt, Config, ConfigList, ConfigConfigList, ConfigConfigConfigList

using Distributions

const ConfigInt = Int16       # Change to UInt8 if N > 127
const Config = Array{ConfigInt,1}   # A config
const ConfigList = Array{Array{ConfigInt,1},1}  # A list of configs correponding to N, K
const ConfigConfigList = Array{Array{Array{ConfigInt,1},1},1}  # A list of lists of configs corresponding to N
const ConfigConfigConfigList = Array{Array{Array{Array{ConfigInt,1},1},1},1} 

const Population = Array{Int64,1}
const PopList = Array{Population,1}

