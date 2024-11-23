export SpacOutput, SpacOutputs


@with_kw mutable struct SpacOutput{FT}
  ## Original Output
  ET::FT = 0
  Tr::FT = 0
  Es::FT = 0
  Ei::FT = 0
  Esb::FT = 0
  RS::FT = 0
  GW::FT = 0
  SM::Vector{FT} = zeros(3)
end


@with_kw mutable struct SpacOutputs{FT}
  ## Original Output
  ntime::Int = 10
  ET::Vector{FT} = zeros(ntime)
  Tr::Vector{FT} = zeros(ntime)
  Es::Vector{FT} = zeros(ntime)
  Ei::Vector{FT} = zeros(ntime)
  Esb::Vector{FT} = zeros(ntime)
  RS::Vector{FT} = zeros(ntime)
  GW::Vector{FT} = zeros(ntime)

  SM::Matrix{FT} = zeros(ntime, 3)
end


function Base.setindex!(res::SpacOutputs, r::SpacOutput, i::Int64)
  fields = fieldnames(SpacOutput)
  @inbounds for f in fields[1:end-1]
    getfield(res, f)[i] = getfield(r, f)
  end
  res.SM[i, :] .= r.SM
  return res
end
