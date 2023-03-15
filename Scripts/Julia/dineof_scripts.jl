function nanmean(x,dim)
    m = isnan.(x)
    x2 = copy.(x)
    x2[m] .= 0
    return sum(x2,dims=dim) ./ sum(.!m,dims=dim)
end

include("dineof_cvp.jl")
include("coverage.jl")
