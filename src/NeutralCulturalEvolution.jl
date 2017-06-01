module NeutralCulturalEvolution
include("../src/aliases.jl")
#include("../src/innovation.jl")
include("../src/innovation_collection.jl")
include("../src/neutral_poplist.jl")
include("../src/conformist_poplist.jl")
include("../src/turnover.jl")
include("../src/nearly_neutral_poplist.jl")
include("../src/freq_scaled_fitness.jl")
include("../src/slatkin_C.jl")
include("../src/watterson.jl")
include("../src/infsites.jl")
include("../ewens/ewens.jl")
include("../ewens/slatkin_J.jl")
include("../ewens/stewart.jl")
include("../ewens/slatpart.jl")
include("../ewens/stirling1.jl")
include("../experiments/dataframe_io.jl")
include("../Rcall/poweRlaw.jl")   # commented out on 12/21/16 since doesn't load on pardosa

end

using NeutralCulturalEvolution
