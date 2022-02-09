
struct Tstruct
		rscale::Vector{Float64}
		Pk::Matrix{Float64}
		Hkc::Matrix{Float64}
		Hks::Matrix{Float64}
end

struct Qstruct
		q::Int64
		Q::Matrix{Float64}
end

struct Rstruct
		cn::Int64
		def_j::Vector{Int64}
		def_col::Matrix{Float64}
		def_R::Matrix{Float64}
		R::Matrix{Float64}
		order::Vector{Int64}
		def_ord::Vector{Int64}
		def_pow::Vector{Int64}
		def_num::Vector{Int64}
end

"""
 The purpose of this function is to compute tilde{R} that is used
 for evaluating the basis functions Psi used in the RBF-QR method.
--- ep (scalar) : The shape parameter
--- xk[1:N,1:2] : The center points in polar coordinates [r,theta] and
                  scaled to the unit disc.
"""
function init_psi_2d(ep::Float64, xk::Matrix{Float64}, rr::Float64, cc::Matrix{Float64})
	# these tolerances can be tuned to improve performance  
  mp = eps() # machine precision
  tolD = 1e4 * mp # The smallest number considered safe to invert
  tol = 2; # if a pivot drops more than this, it is considered suspect
	N = size(xk,1)
	# find the polynomial degree that N basis functions correspond to in n
	# dimensions
	jN = degree(N,2)
	# compute the first part of C, up to degree jN 
  Psi, C, T = addcblocks(ep, xk, jN)	
	Q, R = IncQR(C, ep^2, Psi.j, tol)
	iterate = true
	q = Q.q
	posN = -1
	jmax = jN
	iter = 1
	while iterate
		#--- Add one block to C
    jmax = jmax+1
    Psi, C, T = addcblocks(ep,xk,jmax,Psi,T)
		    #--- If we have less than N columns, then we cannot be done
    if q > N
      if posN < 0 # First time, find out location to compare with
				jN = Psi.j[R.order[N]]
				pos = findall(x -> x .== jN, Psi.j)
				posN = pos[1] # The location of d_j,0
      end		

      #--- Check the magnitude of the latest block. 
      jtest = Psi.j[end]
      pos = findall(x-> x .== jtest, Psi.j)
      pos = pos[1]
      p1 = [1,posN] # First and Nth, either one may be smallest
      p2 = pos
      relsc = EvalD_2D(ep,p1,p2,Psi.j,Psi.m,Psi.p,tolD)
			relsc = maximum(relsc)
	    #--- Block small enough means we are done.
      if relsc * exp(0.223 * jtest + 0.212 - 0.657 * mod(jtest, 2)) < mp
				break
      end
    end
    #--- Update the QR-factorization with the new block
    pos = findall(x -> x == Psi.j[end], Psi.j)
    Q, R = IncQR(C, ep^2, Psi.j[pos], tol, Q, R)
    q = Q.q
		iter = iter + 1
  end
  M = length(R.order)
  Rt = R.R[1:N, 1:N] \ R.R[1:N, N+1:M]
  p1 = R.order[1:N]
	p2 = R.order[(N+1):M]
  if M > N # If there is a part to scale
    D = EvalD_2D(ep, p1, p2, Psi.j, Psi.m, Psi.p, tolD)
    Rt = D .* Rt
  end
	Psi.ep = ep
	Psi.xk = xk
	Psi.Rt = Rt
	Psi.columns = R.order
	Psi.rr = rr
	Psi.cc = cc
	return Psi
end

"""
Compute new blocks of the coefficient matrix C up to degree jmax
--- ep (scalar) : The shape parameter
--- xk[1:N,1:2] : The center points as in InitPsi_2D
--- jmax (scalar) : The degree to stop at
--- Psi (struct) : Supplied the second or higher call
--- T (struct) : Supplied the second or higher call
"""
function addcblocks(ep::Float64, xk::Matrix{Float64}, jmax::Int64, varargs...)
	if length(varargs) == 0
		Psij = zeros(Int64,0)
		Psim = zeros(Int64,0)
		Psip = zeros(Int64,0)
		Psics = zeros(Int64,0)
		j0 = 0
		Pk = ones(size(xk,1), jmax + 1)
		rscale = exp.(-ep^2 .* xk[:,1].^2)
		Hkc = zeros(length(xk[:,2]),0)
		Hks = zeros(length(xk[:,2]),0)
	else
		Psi = varargs[1]
		Psij = Psi.j
		Psim = Psi.m
		Psip = Psi.p
		Psics = Psi.cs
		j0 = Psij[end] + 1
		T = varargs[2]
		Pk = hcat(T.Pk, ones(size(xk,1), jmax -j0 + 1))
		rscale = T.rscale
		Hkc = T.Hkc
		Hks = T.Hks
	end
  j = zeros(Int, 0); m = zeros(Int, 0); p = zeros(Int, 0); odd = mod(j0+1,2)
	@inbounds for k = j0:jmax
		odd = abs(odd - 1)
    append!(p, odd * ones(k + 1))
		append!(j, k * ones(k + 1))
		q = zeros(k + 1)
    q[1:2:k+1] = (0:(k - odd)/2)
    q[2:2:k+1] = abs(odd-1):((k-odd)/2)
    append!(m,q[1:k+1])
	end
	#--- Fill in trigs and powers that will be reused later
	Hkc = hcat(Hkc,cos.(xk[:,2] * collect(max(1,j0):jmax)'))
	Hks = hcat(Hks,sin.(xk[:,2] * collect(max(1,j0):jmax)'))
  @inbounds for k = max(1,j0):jmax
    Pk[:,k+1] = xk[:,1] .* Pk[:,k]
  end
	cs = zeros(Int, size(j)) # find which positions are sine and cosine
	pos = findall(x -> x > 0, 2*m .+ p) # cos = 1, sin = -1
	cs[pos[1:2:end]] .= 1
	cs[pos[2:2:end]] .= -1	

#--- compute the new blocks of the coefficient matrix
  M = length(j) 
  cscale = 2*ones(M)            # Column scaling of C = b_[j,m]
  pos = findall(x -> x == 0, 2*m .+ p);   cscale[pos] .= 0.5 * cscale[pos]
  pos = findall(x -> x == 0, j .- 2*m);   cscale[pos] .= 0.5 * cscale[pos]
  C = Pk[:, j .+ 1]; # the powers of r_k and then the trig part
  pos = findall(x -> x == 1, cs);  C[:,pos] .= C[:,pos] .* Hkc[:,2*m[pos] .+ p[pos]]
  pos = findall(x -> x == -1, cs); C[:,pos] .= C[:,pos] .* Hks[:,2*m[pos] .+ p[pos]]
  C = C .* (rscale * cscale')
  a = (j .- 2*m .+ p .+ 1)./2; b = hcat(j .- 2*m .+ 1, (j .+ 2*m .+ p .+ 2) ./2)
  z = ep^4 .* xk[:,1].^2
  @inbounds for k = 1:M
    C[:,k] = C[:,k] .* hypergeom12(a[k], b[k,:], z)
  end
	Psi =	Psistruct(
	append!(Psij, j),
	append!(Psim, m),
	append!(Psip, p),
	append!(Psics, cs),
	ep,
	xk,
	zeros(0,2),
	zeros(0),
	0,
	zeros(0,2),
	zeros(0,2)
	)		
	T = Tstruct(
		rscale,
		Pk,
		Hkc,
		Hks
	)
	return Psi, C, T
end

function hypergeom12(a::Float64, b::Vector{Float64}, x::Vector{Float64})
#
# In our case, I will only implement the 1F2 case, which we have.
#		
  if (length(a) != 1 || length(b) != 2)
    @error "Wrong number of arguments in hypergeom for 1F2"
  end
  #
  # the first coefficient is 1 always
  # 
  alpha = 1
  pos = findmax(abs.(x)) # Could possibly be complex
  v = ones(size(x))
  y = v
  n = 0
  test = 1;
  mp = eps(); # machine precision
  while (abs(test) > mp)
    #
    # the hypergeometric coefficient is computed from the arguments.
    # http://en.wikipedia.org/wiki/Hypergeometric_function#The_series_pFq
    #
    alpha = alpha*(a+n) / (b[1]+n) / (b[2]+n)
    n += 1
    #
    # power and factorial can also be done recursively.
    #
    v = (1/n) * v.* x
    y = y .+ alpha*v
    test = alpha * v[pos[2]]
	end
	return y
end

"""
 This function performs an incremental QR-factorization
 with a block pivoting strategy specially designed for the RBF-QR
 method.
  
 NOTE: This code is a bit of a hack and needs more work to be robust.

--- newC[1:N,:] : The new columns to incorporate into the QR-fac
--- bf         : The scaling factor between two blocks  
--- bj[:]     : The block number of each column. Columns should come in
                 whole blocks.
--- tol        : Determines when to cut out columns   
--- Q (struct) : Information about Q and pivoting
--- R (struct) : The resulting columns in R
"""
function IncQR(newC::Matrix{Float64}, bf::Float64, bj::Vector{Int64}, tol::Int64, varargs...)
	mp = eps() #Machine precision
	zz = 10*mp
	N = size(newC,1)
# Check if first time operating on C
	if length(varargs) == 0 # First time
		Rcn = 0 # The number of handled/seen columns
		Rdef_j = zeros(0)
		Rdef_col = zeros(N,0)
		Rdef_R = zeros(N,0)
		RR = zeros(N,0)
		Rorder = zeros(0)
		Rdef_ord = zeros(0)
		Rdef_pow = zeros(0)
		Rdef_num = zeros(bj[end])
		Qq = 1 # The next q-vector
		QQ = zeros(N,0)
	elseif length(varargs) == 2
		Q = varargs[1]
		R = varargs[2]
		Rcn = R.cn # The number of handled/seen columns
		Rdef_j = R.def_j
		Rdef_col = R.def_col
		Rdef_R = R.def_R
		RR = R.R
		Rorder = R.order
		Rdef_ord = R.def_ord
		Rdef_pow = R.def_pow
		Rdef_num = R.def_num
		Qq = Q.q # The next q-vector
		QQ = Q.Q
	else
			@error "wrong numer of input arguments"
	end
	cn = size(newC,2)
	cols = Rcn .+ collect(1:cn) # Original column numbers
	Rcn = Rcn + cn     # New starting point for next time
	# --- Compute upper R-values for all new columns
	R0 = QQ' * newC

	# --- If we have not yet filled Q 
	if Qq <= N 
		# --- Remove the Q-component from all columns
		newC = newC .- QQ * (QQ' * newC)
		# --- Work block by block through the matrix (assuming ordered)
		@inbounds for j = bj[1] : bj[end]
			 if (max(1,j) == length(Rdef_num)+1)
				  push!(Rdef_num,0)
			 else
			Rdef_num[max(1,j)] = 0 # Never defer for j=0
				  end
			# --- Locate the current block
			pos = findall(x -> x == j, bj)
			n = length(pos)
			# --- Check for deferred columns
			dpos = findall(x -> x == j, Rdef_j)
			nb = n + length(dpos)
			# --- Store actual column numbers for later use
			order = vcat(cols[pos], Rdef_ord[dpos])
	    # --- QR-factorize the block with pivoting
		  F = qr(hcat(newC[:,pos], Rdef_col[:,dpos]), Val(true)) # A[:,Enew] == Q*R
			Qnew = F.Q
			Rnew = F.R
			Enew = vec(F.p)
			# --- Check the newly computed diagonal elements for significance
		  dR = abs.(diag(Rnew))
	    #RorigR{j+1}=dR;
			# --- Look for elements that drop significantly in magnitude
			mag = log10.(dR)
			ll = length(mag)
			if (ll < nb)
				 append!(mag, mag[ll] * ones(nb))
			end
			diff = mag[1:end-1] .- mag[2:end]
			pp = findall(x -> x > tol, diff) # As default tolerance, I suggest 2 (100 times)
			pp = vcat(pp,nb) # Later used as intervals, adding end 
			# --- Collect the part of R corresponding to old q-vectors
			upper = hcat(R0[:,pos], Rdef_R[1:Qq-1,dpos])

	    # --- Split into two parts. The one to keep and the one to defer  
			upper = upper[:,Enew[1:pp[1]]] # Sort and cut
  
			# ---  Remove the parts that will not be used
			Qnew = Qnew[:, 1:min(N,pp[1])]
			Rnew = Rnew[1: min(N,pp[1]), collect(1:pp[1])]
			if Qq-1+pp[1] > N # We have more columns than we need left
				rows = N-Qq+1
				Rnew = Rnew[1:rows,:]
				Qnew = Qnew[:,1:rows]
			end
			ddR = diag(Rnew);
			if size(ddR,1) == size(ddR,2)
				ddR = ddR[1] # Special if just one row
			end
			#  R.selectedR{j+1} = abs(ddR); 
		  # --- Reorthogonalise the new basis vectors (otherwise drift kills accuracy)
			Qnew = Qnew - QQ * (QQ' * Qnew)
			#--- Renormalize also?
			qnorm = sqrt.(diag(Qnew' * Qnew))
			Qnew = Qnew * diagm(1 ./ qnorm)
			# --- Update Q and R
			RR = hcat(RR, vcat(upper, Rnew, zeros(N-Qq+1-size(Rnew,1),pp[1])))
			QQ = hcat(QQ, Qnew)
			Qq = Qq + size(Qnew,2)
			Rorder = vcat(Rorder, order[Enew[1:pp[1]]])
			 # --- Apply the new q-vectors to the remaining vectors
#			println( Qnew'*newC[:,pos[end]+1:end])
			R0 = vcat(R0, hcat(zeros(size(Qnew, 2), pos[end]), Qnew' * newC[:,pos[end]+1:end]))
			newC[:,pos[end]+1:end] .= newC[:,pos[end]+1:end] - Qnew * (Qnew' * newC[:,pos[end]+1:end]);
			# --- Handle the deferred columns. Assuming tol corresponds to blockdiff
			if length(pp) > 1       # There are columns to defer
				# --- diff = -log10(bf^p) = - p log10(bf)
				if bf < 1
					p = -diff[pp[1:end-1]] ./ log10.(bf) # Default bf = ep^2
				else
					p = ones(size(diff[pp[1:end-1]])) # We can move anything if ep is large
				end
				pow = zeros(pp[end])
				@inbounds for ll = 1:(length(pp)-1)
					pow[(pp[ll]+1):pp[ll+1]] .= floor(sum(p[1:ll]))
				end
				# --- First handle all new deferred columns
				p1 = findall(x -> x .<= n, Enew[pp[1]+1:end])
				p1 = p1 .+ pp[1]
				loc = Enew[p1]
				ll = length(loc)
				Rdef_num[j] = Rdef_num[j] + ll # Counter
				Rdef_col = hcat(Rdef_col, newC[:, pos[loc]])
				Rdef_R = hcat(Rdef_R, vcat(R0[:, pos[loc]], zeros(N-Qq+1, ll)))
				Rdef_ord = vcat(Rdef_ord, cols[pos[loc]])
				Rdef_pow = vcat(Rdef_pow, pow[p1])
				if sum(mag[p1] .<= log10(10*zz)) == length(p1)
#         'Zero'
					pow[p1] .= pow[p1] .+ 5 # Zero, move forward
				end
				if Qq > N
#         'Already at the end'
					pow[p1] .= 1 # Take these into the final R at the end
				end
				Rdef_j = vcat(Rdef_j, j .+ pow[p1])
				# --- Go through old columns to redefer
				p2 = findall(x -> x .> n, Enew[pp[1]+1:end])
				p2 = p2 .+ pp[1]
				@inbounds for k = 1:length(p2) 
					loc = Enew[p2[k]] - n
					# --- Check if we can move it even further. It is small also here
#         'Redefer?'
					pow[p2[k]] = pow[p2[k]] - Rdef_pow[dpos[loc]]
					if mag[p2[k]] <= log10(10*zz)
#           'Zero'
						pow[p2[k]] = 5 # Zero, move forward
					end
					if pow[p2[k]] > 0
#           'Moving'
						Rdef_pow[dpos[loc]] = Rdef_pow[dpos[loc]] + pow[p2[k]]
						Rdef_j[dpos[loc]] = j + pow[p2[k]]
					else
    #				 R.def_pow(dpos(loc))
					  powerval = pow[p2[k]]
						magval = mag[p2[k]]
						# MAY WANT TO CHANGE THIS BACK TO AN ERROR
						@warn "Found an undeferable column. This case is not defined yet"
						# Probably we should then decide to add the column
						pow[p2[k]] = 1
						Rdef_pow[dpos[loc]] = Rdef_pow[dpos[loc]] + pow[p2[k]]
						Rdef_j[dpos[loc]] = j + pow[p2[k]]
					end        
				end
			end
			# --- Apply the new Q to the deferred columns
			nq = size(Qnew,2)
			Rdef_R[Qq.-nq:Qq-1,:] = Qnew' * Rdef_col
			Rdef_col = Rdef_col .- Qnew * (Qnew' * Rdef_col)
	    if size(RR,2) >= N
		    # --- Put the deferred columns back into R at the end.
			  pos = findall(x -> x .> j, Rdef_j)
			 # --- Make nearly zero into exactly zero
				@inbounds for q = 1:length(pos)
					where_idx = findall(x -> x .<= 10*zz, abs.(Rdef_R[:,pos[q]]))
					Rdef_R[where_idx,pos[q]] .= 0
				end
				RR = hcat(RR, Rdef_R[:,pos]) # possibly sort them
				Rorder = vcat(Rorder, Rdef_ord[pos])
				# --- Put in the columns that were not treated yet.
				pos = findall(x -> x .> j, bj)
				rz = size(R0,1)
				RR = hcat(RR, vcat(R0[:,pos], (QQ[:,rz+1:end]' * newC[:,pos])))
				Rorder = vcat(Rorder, cols[pos])
				break
			end
		end
	else
		# --- If we have filled all columns of Q, just project and join
		RR = hcat(RR, R0)
		Rorder = vcat(Rorder, cols)
	end
	#origR = RorigR
	#selectedR = RselectedR
	Q = Qstruct(
				Qq,
				QQ
			)
	R = Rstruct(
				Int64.(Rcn),
				Rdef_j,
				Rdef_col,
				Rdef_R,
				RR,
				Rorder,
				Rdef_ord,
				Rdef_pow,
				Rdef_num
			)
	return Q, R
end

"""
The purpose of this function is to compute the scaling effect of 
 D_1^{-1} and D_2 applied to the correction matrix in the RBF-QR method.
--- ep (scalar) : The shape parameter
--- p1 (vector) : Indices for elements in D_1 
--- p2 (vector) : Indices for elements in D_2
--- j,m,p (vectors) : Identifiers for the expansion functions T_{j,m}  
"""
function EvalD_2D(ep,p1,p2,j,m,p,tol)
  tovec(f) = typeof(f) == Vector{Int64} ? f : [f]
  D = zeros(length(p1),length(p2))
  ep2 = 0.5*ep*ep
  pmin = maximum(j[p1]); # Largest negative power that could occur
  if (ep2 < 1) # May need to limit the powers
    if (ep2 > tol)
      pmin = min(pmin, Int64.(floor(log(tol)/log(ep2))))
    else    
      pmin = 0
    end
  end
#--- Precompute all positive and negative powers of ep that are present.
  pmax = maximum(j[p2]);
	epp = zeros(1 + pmin + pmax)
  epp[0+pmin+1] = 1
  @inbounds for pp = 1:pmax
    epp[pp+pmin+1] = ep2*epp[pp + pmin]
  end
  ep2 = 1/ep2
  @inbounds for pp = 1:pmin
    epp[-pp+pmin+1] = ep2*epp[-pp+pmin+2]
  end  
#--- Precompute powers of 2, both negative and postive
  mmax = maximum(m[p2])
  mmin = maximum(m[p1]) # Largest negative
	twop = zeros(1 + mmin + mmax)
  twop[0 + mmin + 1] = 1  # 2^0
  @inbounds for pp = 1:mmax
    twop[pp + mmin + 1] = 4*twop[pp+mmin]
  end
  @inbounds for pp = 1:mmin
    twop[-pp+mmin+1] = 0.25*twop[-pp+mmin+2] 
  end
#--- The column values stay the same for each row
  powj = j[p2]
  powm = m[p2]
  @inbounds for k = 1:length(p1) # For each row
    pow = tovec(powj .- j[p1[k]])
#--- Remove too negative powers
    pos = findall(x -> x >= -pmin, pow)
		D[k,pos] .= epp[pow[pos] .+ pmin .+ 1]
    pow = powm .- m[p1[k]]
    D[k,:] .= D[k,:] .* twop[pow .+ mmin .+ 1]
  end
#--- This part does the ratios of factorials
  f1 = tovec(Int64.((j[p2] .+ 2 .* m[p2] .+ p[p2]) ./ 2))
  f2 = tovec(Int64.((j[p2] .- 2 .* m[p2] .- p[p2]) ./ 2))
  f3 = tovec(Int64.((j[p1] .+ 2 .* m[p1] .+ p[p1]) ./ 2))
  f4 = tovec(Int64.((j[p1] .- 2 .* m[p1] .- p[p1]) ./ 2))
  fmax = max(maximum(f1), maximum(f3))
  fp = cumprod(vcat(1, 1:fmax)) # Because 0!=1

  @inbounds for k = 1:length(p1) # For each row
    numer = ones(1, size(D,2))
    denom =  ones(1, size(D,2))
    pos = findall(x -> x .> f3[k], f1)
		if !isempty(pos)
	    denom[pos] .= (1 / fp[f3[k] + 1]) * fp[f1[pos] .+ 1]
		end
    pos = findall(x -> x .> f4[k], f2)
		if !isempty(pos)
			denom[pos] .= (1 / fp[f4[k] + 1]) * denom[pos] .* fp[f2[pos] .+ 1]
		end
		pos = findall(x -> x .< f3[k], f1)
		if !isempty(pos)
			numer[pos] .= fp[f3[k].+1] ./ fp[f1[pos].+1]
		end
    pos = findall(x -> x .< f4[k], f2)
		if !isempty(pos)
	    numer[pos] .= numer[pos] .* (fp[f4[k] .+ 1] ./ fp[f2[pos] .+ 1])
		end
    D[k,:] .= D[k,:] .* numer[:] ./ denom[:]
  end
	return D

end

"""
Find the polynomial degree that N basis functions correspond to in n
dimensions
"""
function degree(N::Int64, n::Int64)
	# dim: Find the dimension of the polynomial space of degree K in n dimensions
	dim(K,n) = prod(((K+1):(K+n))./(1:n))
	K = 0.0
  @inbounds for k = 0:N-1 # K(N) cannot be larger than N-1 (1D-case)
    if (dim(k,n) >= N)
      K = k
      break
    end
  end
	return K
end
