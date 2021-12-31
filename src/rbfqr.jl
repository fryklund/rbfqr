module rbfqr

using LinearAlgebra

mutable struct Psistruct
	j::Vector{Int64}
	m::Vector{Int64}
  p::Vector{Int64}
	cs::Vector{Int64}
  ep::Float64
	xk::Matrix{Float64}        
  Rt::Matrix{Float64}
	columns::Vector{Float64}
	rr::Float64
	cc::Matrix{Float64}
	A0::Matrix{Float64}
end

include("RBFdiscretization.jl")
include("rbfqr_diffmat_2d.jl")
include("init_psi_2d.jl")
include("../tests/rbfqr2dTests.jl")

export rbfqr_diffmat_2d
export Psistruct


end # end


