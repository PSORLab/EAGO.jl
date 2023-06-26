# Added from Oscar Dowson's code for JuMP at https://github.com/jump-dev/JuMP.jl/pull/2484
generate_precompile = true
package_path = ""
example_path = ""

module Foo

    using JuMP
    using EAGO

    function stress_precompile()
        for file in readdir("precompiles")
            if !endswith(file, ".jl")
                continue
            end
            include(file)
        end
        return
    end

    if generate_precompile
        using SnoopCompile
        tinf = @snoopi_deep Foo.stress_precompile()
        ttot, pcs = SnoopCompile.parcel(tinf)
        SnoopCompile.write("precompiles", pcs)
        for file in readdir("precompiles")
            if !endswith(file, ".jl")
                continue
            end
            src = joinpath("precompiles", file)
            m = match(r"precompile\_(.+)\.jl", file)
            modules = split(m[1], ".")
            modules = vcat(modules[1], "src", modules[2:end])
            if !(modules[1] in ["EAGO"])
                continue
            end
            dest = joinpath(package_path, modules..., "precompile.jl")
            @show dest
            cp(src, dest; force = true)
        end
    end
end