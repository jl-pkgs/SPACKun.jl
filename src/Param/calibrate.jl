ParRange = [
  (; name=:rs, x0=[70, 100, 200.0], lower=[10.0, 40.0, 40.0], upper=[1000, 1000, 1000.0]), # 气孔导度阻力
  (; name=:Hc, x0=[0.12, 0.5, 1.0], lower=[0.05, 0.1, 0.1], upper=[0.5, 1.0, 2.5])         # 冠层高度
]

# get field of param, and combine into a vector
unlist_field(l::Vector, key::Symbol) = vcat(map(x -> getfield(x, key), l)...)

select_param(parNames) = filter(x -> x.name ∈ parNames, ParRange)

# set parameters in soil
function set_param!(param, theta, par)
  k = 0
  for i = eachindex(par)
    d = par[i]
    n = length(d.x0)
    inds = k+1:k+n
    p = getfield(param, d.name)
    p .= theta[inds]
    k += n
  end
end

function init_param(::Nothing=nothing; soiltype=2, lc=11, perc_full=0.7, kw...)
  soil = Soil{Float64}(; soiltype, lc)
  soil.θ .= soil.param.θ_sat * perc_full # 70% of θ_sat
  soil
end

function init_param(theta::Vector{Float64}; par=ParRange, soiltype=2, lc=11, kw...)
  soil = init_param(nothing; soiltype, lc, kw...)
  set_param!(soil.param, theta, par)
  soil
end


export ParRange, unlist_field, select_param, set_param!, init_param
