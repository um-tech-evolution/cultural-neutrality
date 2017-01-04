module NeutralCulturalEvolution
include("aliases.jl")
include("neutral_poplist.jl")
include("conformist_poplist.jl")
include("nearly_neutral_poplist.jl")
include("freq_scaled_fitness.jl")
include("slatkin_C.jl")
include("watterson.jl")
include("../ewens/ewens.jl")
include("../ewens/slatkin_J.jl")
include("../ewens/stewart.jl")
include("../ewens/slatpart.jl")
include("../ewens/stirling1.jl")
include("../Rcall/poweRlaw.jl")   # commented out on 12/21/16 since doesn't load on pardosa

end

using NeutralCulturalEvolution
