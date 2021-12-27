module RBFdiscretization 

using LinearAlgebra

export ndgrid
export halton


function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = (reshape(collect(v1), 1, m))
    v2 = (reshape(collect(v2), n, 1))
    return (repeat(v1, n, 1), repeat(v2, 1, m))
end		


"""
This routine creates the identical sequences as the routine
haltonseq that is available from Matlab Central
numpts  (scalar) number of points to generate in Halton sequence
ndims   (scalar) number of dimensions, should be <=6
"""
function halton(numpts::Int64, ndims::Int64)
	if ndims > 6
			@error "ndims > 6"
	end
	p = [2,3,5,7,11,13]
	H = zeros(numpts, ndims)
	for k = 1:ndims
			N = p[k]; v1 = 0.0; v2 = 0:N-1; lv1 = 1;
			while lv1 <= numpts
					v2 = v2[1:max(2, min(N,Int64(ceil((numpts+1)/lv1))))] ./ N
					x1, x2 = ndgrid(collect(v2),vec(collect(v1)))
					v1 = x1 .+ x2; v1 = v1[:]; lv1 = length(v1)
			end
			H[:,k] = v1[2:numpts+1]
	end
	return H
end

end # module
