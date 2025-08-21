module Hadron
include("bootstrappReturnType.jl")
include( "bootstrapp.jl")

include("io.jl")
include("visualisation.jl")
export Bootstrap, ts_boot, boot, describeBoot, show

export histboot,analyse,histboot!

# Write your package code here.
function __init__()
    @static if !isdefined(Base, :get_extension)
        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
            include("../ext/MakieExtension.jl")
        end
    end

    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        if exc.f in [histboot,analyse]
            if isempty(methods(exc.f))
                print(io, "\n$(exc.f) has no methods, yet. Makie has to be loaded for the plotting extension to be activated. Run `using Makie`, `using CairoMakie`, `using GLMakie` or any other package that also loads Makie.")
            end
        end
    end
end
end
