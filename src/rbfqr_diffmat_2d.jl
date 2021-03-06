
struct Pstruct
	    F::Matrix{Float64}
      Q::Matrix{Float64}
			G::Matrix{Float64}
			S::Matrix{Float64}
			R::Matrix{Float64}
		d2Te::Matrix{Float64}
    deg::Int64
     Te::Matrix{Float64}
    Hec::Matrix{Float64}
    Hes::Matrix{Float64}
     xe::Matrix{Float64}
    re2::Vector{Float64}
     Pe::Matrix{Float64}
    rsc::Vector{Float64}
      r::Vector{Float64}
end



"""
 [A,Psi]=RBF_QR_diffmat(op,xe,xk,ep) % First time 
 [A,Psi]=RBF_QR_diffmat(op,xe,Psi)   % If Psi is already computed

 Computes a differentiation matrix (weights for RBF-FD stencils)
 using Gaussian RBFs at evaluation point(s) xe for RBFs placed at
 the nodes xk. 

--- op (char)   : Alternatives '1', 'x', 'y', 'xx', 'xy', 'yy' ,'L',
                   'Ln' where n is a number indicating the degree
                   of the Laplacian.
--- xe(1:M,1:2)  : The evaluation point(s) in Cartesian coordinates
--- xk(1:N,1:2)  : The node points (no particular scaling assumed)
--- ep (scalar)  : The (constant) shape parameter 
--- Psi (struct) : Generated by a call to this function. Defines
                   the RBF-QR basisfunctions Psi.

--- Check if Psi exists. Otherwise, scale nodes and compute Psi.
"""
function rbfqr_diffmat_2d(op, xe::Matrix{Float64}, varargs...)
# --- Check if Psi exists. Otherwise, scale nodes and compute Psi.
	if length(varargs) >= 2
		xk = varargs[1]
		ep = varargs[2]
	  # If only one evaluation point, place that in zero.
		if (size(xe,1)==1)
			cc = xe
		else  
			cc = sum(xk,dims=1) ./ size(xk,1) # Center of disc
		end  

		xk = hcat(xk[:,1] .- cc[1], xk[:,2] .- cc[2])
		xe = hcat(xe[:,1] .- cc[1], xe[:,2] .- cc[2])
		r = sqrt.(sum(xk.^2, dims = 2))
		re = sqrt.(sum(xe.^2, dims = 2))
		rr = max(maximum(r),maximum(re))          # Radius of disk
    
		xk = hcat((1/rr) .* r, atan.(xk[:,2], xk[:,1])) # Polar coordinates
		xe = hcat((1/rr) .* re, atan.(xe[:,2] ,xe[:,1])) # Polar coordinates
		ep = ep*rr
		Psi = init_psi_2d(ep, xk, rr, cc)
		#--- Also compute the interpolation matrix which is reused for all ops 
		A0, P = rbf_qr_mat_2d(Psi, "1", Psi.xk)
		Psi.A0 =	A0
	elseif length(varargs) == 1
		Psi = varargs[1]
		xe = hcat(xe[:,1] .- Psi.cc[1], xe[:,2] .- Psi.cc[2])
		re = sqrt.(sum(xe.^2, dims=2))
		xe = hcat((1 / Psi.rr) .* re, atan.(xe[:,2], xe[:,1])); # Polar coordinates  
	end
	# --- Compute the differentiation matrix/ces
	if typeof(op) != Vector{String}
	  #--- We are just computing one operator
		A, P = rbf_qr_mat_2d(Psi, op, xe)
		A = A	/ Psi.A0
		A = rescale_op(A, Psi.rr, op)
	else
		numop = length(op)
		A = Array{Any}(undef,numop)
			# TODO: dictionary
		@inbounds for i = 1:numop
			if i == 1
				var = xe
			else
				var = P
			end
			println(typeof(A[i]))
			A[i], P .= rbf_qr_mat_2d(Psi,op[i],var)
			A[i] .= A[i] / Psi.A0
			A[i] .= rescale_op(A[i], Psi.rr, op[i])
		end
	end
	return A, Psi
end

"""
 The purpose of this function is to return an evaluation matrix A
 for the given operator op. The elements A_{ij} are given by
 op(psi_j[xe[i]]), where psi is the RBF-QR basis and xe are the
 evaluation points.
--- Psi (struct) : Defines the basis functions
--- op (string) : Defines the operator
--- var (double(1:N,1:2) or struct) : Either the evaluation points
---    xe or a structure containing precomputed data including xe

--- Call patterns
--- A, P = rbf_qr_mat_2d(Psi,op,xe) First time for these points (xe)
--- A, P = rbf_qr_mat_2d(Psi,op,P)  Subsequent times with same xe
"""
function rbf_qr_mat_2d(Psi::Psistruct, op::String, varargs...)
#--- Find out which call sign we got
	if typeof(varargs[1]) == Pstruct
		P = varargs[1]
		xe = P.xe
	else
		P = nothing
		xe = varargs[1]
	end
	deg, diff, op = RBF_QR_parse(op)
#--- Precompute functions that are needed for evaluation of T_{j,m}
  ep = Psi.ep
  Rt = Psi.Rt
  N = size(Psi.Rt,1)
  order = Int64.(Psi.columns)
  M = length(order)
  j = Psi.j[order]
  m = Psi.m[order]
  cs = Psi.cs[order]
  p = Psi.p[order]
  P = RBF_QR_precomp_2D(j,m,p,ep,xe,deg,P)
	#--- Evaluate the basis functions
	if deg == 0
    V = (P.rsc*ones(1,M)) .* P.Pe[:, m .+ 1] .* P.Te[:, j .- 2*m .+ 1]
    pos = findall(x -> x == 1, cs)
		V[:,pos] .= V[:,pos] .* P.Hec[:, 2*m[pos] .+ p[pos]]
    pos = findall(x -> x == -1, cs)
		V[:,pos] .= V[:,pos] .* P.Hes[:,2*m[pos] .+ p[pos]]
    A = V[:,1:N] .+ V[:,N+1:end] * transpose(Rt)
	elseif deg == 1
    co1 = hcat(P.Hec[:,1], P.Hes[:,1])
    co2 = [1, -1]
    if op == "x"
      i1 = 1
			i2 = 2
    elseif op == "y"
      i1 = 2
			i2 = 1
    end  
    A = (co1[:,i1] * ones(1,M)) .* P.G
    V = (co1[:,i2] * ones(1,M)) .* P.F
    pos = findall(x -> x == 1, cs)
    A[:,pos] .= A[:,pos] .* P.Hec[:, 2 * m[pos] .+ p[pos]] .+
				co2[i1] * V[:,pos] .* P.Hes[:,2 * m[pos] .+ p[pos]]
    pos = findall(x -> x == -1, cs)
    A[:,pos] .= A[:,pos] .* P.Hes[:,2 * m[pos] .+ p[pos]] .+
				co2[i2] * V[:,pos] .* P.Hec[:, 2 *m[pos] .+ p[pos]]
    
    A = A[:,1:N] + A[:,N+1:end]*transpose(Rt)
	elseif deg == 2
		if op[1] == 'L'
			A = P.Q .+ P.S
			pos = findall(x -> x == 1, cs)
			A[:,pos] .= A[:,pos] .* P.Hec[:, 2 * m[pos] .+ p[pos]]
			pos = findall(x -> x == -1, cs)
			A[:,pos] .= A[:,pos] .* P.Hes[:,2 * m[pos] .+ p[pos]]
		else
			# g(theta)=cos((2m+p)theta) => h(theta)=-sin()
			# g(theta)=sin((2m+p)theta) => h(theta)= cos()
			co1 = hcat(P.Hec[:,1].^2, P.Hes[:,1].^2, P.Hec[:,1] .* P.Hes[:,1], P.Hec[:,2])
			co2 = [1, 1, -1]
			co3 = [2, -2, -1]
			if op == "xx"
				i1 = 1
				i2 = 2
				i3 = 3
			elseif op == "xy"
				i1 = 3
				i2 = 3
				i3 = 4
			elseif op == "yy"
				i1 = 2
				i2 = 1
				i3 = 3
			end  
			A = (co1[:,i1] * ones(1,M)) .* P.Q .+ (co2[i1] * co1[:,i2] * ones(1,M)) .* P.S
			V = (co3[i1] * co1[:,i3] * ones(1,M)) .* P.R

			pos = findall(x -> x == 1, cs)
			A[:,pos] .= A[:,pos] .* P.Hec[:, 2 * m[pos] .+ p[pos]] .-
					V[:,pos] .* P.Hes[:, 2 * m[pos] .+ p[pos]]
			pos = findall(x -> x == -1, cs)
			A[:,pos] .= A[:,pos] .* P.Hes[:,2 * m[pos] .+ p[pos]] .+
					V[:,pos] .* P.Hec[:, 2 * m[pos] .+ p[pos]]
		end	
		A = A[:,1:N] .+ A[:,N+1:end] * transpose(Rt)
	end
	return A, P
end


"""
The purpose of the function is to expand the operator information
 into the needed details
--- op (string[:]) : Describes the operator to evaluate
--- deg (scalar) : The degree of the operator
--- diff (1:dim) : The number of diffs in each dimension
"""
function RBF_QR_parse(op::String)
	#--- Trim away spaces
	pos = findall(x -> x != " ", op)
	op = op[pos]
	deg = Int64[]
	diff = zeros(Int64,3)
	if op[1] == '1'
		deg = 0
	elseif op[1] == 'x' || op[1] == 'y' || op[1] == 'z'
		deg = length(op)
		if deg>2
			error("Mixed derivatives of higher degree than 2 are not implemented")
		end  
		pos = findall(x -> x == 'x', op)
		if length(pos) > 0
			diff[1] = length(pos)
		end
		pos = findall(x -> x == 'y', op)
		if length(pos) > 0
			diff[2] = length(pos)
		end
		pos = findall(x -> x == 'z', op)
		if length(pos) > 0
			diff[3] = length(pos)
		end
	elseif op[1] == 'L'
		diff = Int64[]
		deg = 2
		if length(op) > 1
			deg = 2 * parse(Int64, op[2:end])
		end
		if deg > 10
			error("Hyperviscosity implemented only up to degree L^5")
		end  
	end
	return deg, diff, op
end


function RBF_QR_precomp_2D(j, m, p, ep, xe, deg, P)
	#--- Extract descriptors for the Psi-functions
  tol = 10 * eps() # 10 times the machine precision
  jmax = (j[end])
	pmax = p[end]
  Ne = size(xe,1)
  M = length(j)
  if isnothing(P) # First time
    #--- First compute the basic functions that T are built from
    PF = zeros(Ne, M)
		PQ = zeros(Ne, M) # Used for testing later
		PG = zeros(Ne, M)
		PS = zeros(Ne, M)
		Pd2Te = zeros(0, 2)
		PR = zeros(0,2)
    Pdeg = 0
    PTe = cos.(acos.(xe[:,1]) * (0:jmax)')
			
    PHec = cos.(xe[:,2] * (1:jmax)')
    PHes = sin.(xe[:,2] * (1:jmax)')
  
    Pxe = xe
    Pre2 = xe[:,1] .^ 2  # Only even powers are needed for evaluation points
    PPe = ones(Ne, div(jmax-pmax,2) + 1)
    @inbounds for pp = 1:div(jmax-pmax,2)
      PPe[:,pp+1] = Pre2.*PPe[:,pp]
    end  
    Prsc = exp.(-ep^2 .* Pre2)
    Pr = xe[:,1]
  end
  
  if deg >= 1 && Pdeg < 1
    Pdeg = 1
    Pfac  = 1 .- Pre2     
    pos = findall(x -> x .<= tol, abs.(Pr .- 1))
		Pfac[pos] .= 1
		Pfac = 1 ./ Pfac
    PdTe = .-((Pr .* Pfac) .* (0:jmax)') .* PTe
    PdTe[:,2:end] .= PdTe[:,2:end] .+ (Pfac * (1:jmax)') .* PTe[:,1:end-1]
    PdTe[pos,:] .= ones(length(pos)) * ((0:jmax)').^2
  end

  if  deg >= 2 && Pdeg < 2
    Pdeg = 2
    jj = (0:jmax)'
    pos = findall(x -> x .<= tol, abs.(Pr .- 1)) 

    Pd2Te = .- Pfac .* jj.^2 .* PTe .+ (xe[:,1] .* Pfac) .* ones(1,length(jj)) .* PdTe
    Pd2Te[pos,:] .= ones(length(pos),1)*(jj.^2 .* (jj.^2 .-1)./3)
  end
      
  if deg == 1 && norm(PF) == 0 # Needed now, but not computed yet
    pr0 = findall(x -> x .<= tol, Pr)
    pm0 = findall(x -> x == 0, m) # For s = 0, the coeff is zero
    pm = findall(x -> x .> 0, m)
    #--- (2m+p)T_{j,m}/r=(2m+p)exp(-ep^2r^2)r^{2m-1}T_{j-2m}(r)
    PG[:,pm] .= (Prsc.*Pr)*ones(1,length(pm)) .* PPe[:,m[pm] .+ 1 .- 1].*((-2*ep^2 * Pre2 .* ones(1,length(pm)) .+ 2*ones(Ne,1)*m[pm]') .* PTe[:,j[pm] .- 2*m[pm] .+ 1] .+ Pr*ones(1,length(pm)).*PdTe[:,j[pm] .- 2*m[pm] .+ 1])
			
    PG[:,pm0] .= Prsc*ones(1,length(pm0)).*
			( -2*ep^2*Pr*ones(1,length(pm0)).*PTe[:,j[pm0] .+ 1] .+
				PdTe[:,j[pm0] .+ 1])

    PF[:,pm] .= (Prsc.*Pr) * (2 * m[pm] .+ p[pm])' .* PPe[:, m[pm] .+ 1 .- 1] .* PTe[:, j[pm] .- 2*m[pm] .+ 1]

    PF[:,pm0] .= (Prsc ./ Pr) * ( 2 * m[pm0] .+ p[pm0])' .* PTe[:,j[pm0] .+ 1]

    PF[pr0,pm0] .= ones(length(pr0),1) * (j[pm0] .* cos.((j[pm0] .- 1) ./ 2 * pi))'
  end
  
  if  deg == 2 && norm(PQ) == 0    
    pr0 = findall(x -> x <= tol, Pr)
    pm0 = findall(x -> x ==0, m) # For s=0, the coeff is zero
    pm = findall(x -> x > 0, m)

    cT = 2*ep^2*(2*ep^2*Pre2*ones(1,M) .- 4*ones(Ne,1) * m' .- 1 )
    cT[:,pm] .= cT[:,pm] .* (Pre2 * ones(1,length(pm))) .+ ones(Ne,1) * (2 * m[pm] .* (2 * m[pm] .- 1))'
    %
    cdT = .-4 * ep^2 * Pr*ones(1,M)
    cdT[:,pm] .= cdT[:,pm] .* (Pre2 * ones(1,length(pm))) .+ 4*Pr*m[pm]'
    
    cd2T = ones(Ne,M)
    cd2T[:,pm] .= cd2T[:,pm].*(Pre2*ones(1,length(pm)))
    
    PQ = cT .* PTe[:,j .- 2*m .+ 1] .+ cdT .* PdTe[:,j .- 2*m .+ 1] + cd2T .* Pd2Te[:,j .- 2*m .+ 1]
    PQ .= (Prsc*ones(1,M)) .* PQ
    PQ[:,pm] .= PPe[:, m[pm] .- 1 .+ 1] .* PQ[:,pm]
      
    #--- The rest of the terms have special cases also for m=1
    rinv = Pr
		rinv[pr0] .= 1
		rinv = 1 ./ rinv # Max value is 14.
    PR = zeros(Ne,M)
      
    cT = ones(Ne,1) * (1 .- 2*m)'
    cT[:,pm] .= cT[:,pm] .+ 2 * ep^2 * (Pre2 * ones(1,length(pm)))
    cT[:,pm0] .= cT[:,pm0] .* (rinv.^2 * ones(1,length(pm0))) .+ 2*ep^2
      
    cdT = zeros(Ne,M)
    cdT[:,pm] .= -Pr*ones(1,length(pm))
    cdT[:,pm0] .= .- rinv * ones(1,length(pm0))

    PR = cT .* PTe[:,j .- 2*m .+ 1] .+ cdT .* PdTe[:,j .- 2*m .+ 1]
    PR[pr0,pm0] .= 0; # Limit for m=r=0 is zero in all cases
    PR = (Prsc*(2*m .+ p)') .* PR # This makes the m=p=0 case = 0
    PR[:,pm] .= PR[:,pm] .* PPe[:,m[pm] .- 1 .+ 1]
          
      
    cT = ones(Ne,1)*(2*m.- (2*m .+ p).^2)'
    cT[:,pm] .= cT[:,pm] .- 2*ep^2 * Pre2*ones(1,length(pm))
    cT[:,pm0] .= cT[:,pm0] .* (rinv.^2 * ones(1,length(pm0))) .- 2*ep^2
    
    cdT = zeros(Ne,M)
    cdT[:,pm] .= Pr*ones(1,length(pm))
    cdT[:,pm0] .= rinv*ones(1,length(pm0))
    
    PS = cT .* PTe[:,j .- 2*m .+ 1] .+ cdT .* PdTe[:,j .- 2*m .+ 1]
    PS = (Prsc * ones(1,M)) .* PS
    PS[:,pm] .= PS[:,pm] .* PPe[:,m[pm] .- 1 .+ 1]
    pm00 = intersect(findall(x -> x ==0, m), findall(x -> x ==0, p))
    PS[pr0,pm0] .= 0 # First make all zero
    PS[pr0,pm00] .= ones(length(pr0),1) *
				((-1).^(j[pm00] ./2 .+ 1) .* (2*ep^2 .+ j[pm00].^2))'
	end
		return Pstruct(
	    PF,
      PQ,
			PG,
			PS,
			PR,
	 Pd2Te,
    Pdeg,
     PTe,
    PHec,
    PHes,
     Pxe,
    Pre2,
     PPe,
    Prsc,
      Pr
		)

end



function rescale_op(A,rr,op)
	#--- Adjust the scaling according to derivative
	deg, diff, op = RBF_QR_parse(op)
	return (1/rr).^deg * A
end

function myfindall(f, a::Array{T, N}) where {T, N}
    j = 1
    b = Vector{Int}(undef, length(a))
    @inbounds for i in eachindex(a)
        if f(a[i])
            b[j] = i
            j += 1
        end
    end
    resize!(b, j-1)
    sizehint!(b, length(b))
    return b
end
