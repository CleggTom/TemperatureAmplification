#define parameter structure
struct params
    r::Vector{Float64}
    R::Vector{Float64}
    a::Array{Float64,2}
    Nsp::Int8
end

#define GLV model
function dC!(dC,C,p,t)
    for i = 1:p.Nsp
        dC[i] = C[i] * p.r[i]
        for j = 1:p.Nsp
            dC[i] += C[i] * C[j] * p.a[i,j]
        end
    end
end

#average biomass
#for non-meanfield
function uC_equi(p::params)
    a = [p.a[i,i] for i = 1:p.Nsp]
    ψ = sum(p.a,dims=2) .- a
    -(mean(p.r) / mean(a)) * ( 1 / (1 + ( mean(ψ)/mean(a) ) ) )
end

#thermal sensitivity
function E_R_eco(R,C,ER,EC)
    top = R .* C .* (ER .+ EC)
    bot = R .* C 
    sum(top) / sum(bot)
end