export Population, PopList, ConfigInt, Config, ConfigList, ConfigConfigList, ConfigConfigConfigList

typealias ConfigInt UInt8     # Change to UInt8 if N > 127
typealias Config Array{ConfigInt,1}   # A config
typealias ConfigList Array{Array{ConfigInt,1},1}  # A list of configs correponding to N, K
typealias ConfigConfigList Array{Array{Array{ConfigInt,1},1},1}  # A list of lists of configs corresponding to N
typealias ConfigConfigConfigList Array{Array{Array{Array{ConfigInt,1},1},1},1} 

typealias Population Array{Int64,1}
typealias PopList Array{Population,1}

