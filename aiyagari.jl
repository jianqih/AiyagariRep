begin
    using Parameters # enable @unpack
    using Distributions # enable cdf of Normal
    using Roots         # root finding routines: brent
    using PlutoUI    # printing stuff in the Pluto notebook
    using Plots 
    using Interpolations
end


include("tauchen.jl")
include("interpolation.jl")

include("para.jl")

function coordGridIntp!(xgrid::Array{Float64,1}, xtarget::Array{Float64,1}, 		ibelow::Array{Int64,1}, iweight::Array{Float64,1}; robust = false)
		xg = length(xgrid); #xt = length(xtarget)
		xi  = 1; xlow = xgrid[1]; xhigh = xgrid[2] # initialize
		@inbounds for (it, xtval) in enumerate(xtarget)
			while xi < xg - 1 # find grid point
				if xhigh >= xtval; break; end
				xi += 1
				xlow = xhigh
				xhigh = xgrid[xi + 1]
				end
			iweight[it] = (xhigh - xtval) / (xhigh - xlow) #  careful with last point. if xhigh<xtval might run into problems.
			ibelow[it] = xi
		end
		if robust == true; iweight =  min.(max.(iweight,0.0), 1.0); end
	return iweight, ibelow
end



param = setPar() # check the parameters above

function solveHHproblemEGM(param, r, w)
    @unpack nA, nS,  gA, gS, transS, β, γ, δ = param

    # Define parameters
    maxiter = 2000
    tol = 1E-5 # tolerance
	
	# Useful matrices
    gAA = repeat(gA, 1, nS) # nA x nS matrix of asset grid.
    gSS = repeat(gS', nA, 1) # nA x nS matrix of labor grid. Note the transpose.
    gYY = w.*gSS .+ (1.0+r).*gAA # nA x nS matrix of cash on hand grid.

    # Pre-allocate value and policy funcctions
    cp = gYY - repeat(gA, 1, nS) # consumption policy
    ap = copy(cp)   # asset policy
	endogAA = copy(cp)  # endogenous asset grid 


	## === Function that runs one iteration of EGM
	function endoGridIter(cp)
		expect = β*(1.0+r)* cp.^-γ*transS' # Right hand side of Euler Eq.
		c = expect.^(-1.0/γ)    # Invert marginal util to get contemporaneous C
		endogAA = (c + gAA - w.*gSS)/(1+r) # compute asset state on endogenous grid (note that gAA is the policy function, which is on-grid)

		# Now compute the state on-grid interpolating the endogenous grid
		apNew = zeros(nA, nS)
		for is = 1:nS
			apNew[:, is] = linearGridIntp!(endogAA[:, is], gA, gA, apNew[:, is])

			# Must account for the binding constraint!
			for ia = 1:nA
				if apNew[ia, is] < gA[1] 
					apNew[ia, is] = gA[1]
				else
					break # exploit monotinicity of ia.
				end
			end
		end

		cpNew = gYY - apNew # get updated consumption policy 

		return cpNew, apNew, endogAA
	end
	
	
	## ===== Start iteration ====
	for iter = 1:maxiter
		cpNew, ap, endogAA = endoGridIter(cp)

		# ===== Check if policuufunction has converged =====
        d = maximum(abs.(cp - cpNew))
		cp = cpNew

        if d < tol
			println("Tol. achieved: $d")
			break # break the loop in case we found the Policy Function!
		end

		if iter == maxiter;
        	println("Max iterations achieved. VF did not converge")
        end

		#println("Iter: $iter")
		#println("Tol: $d") # ps these prints will not show in the notebook.
	end

	return (ap = ap,  cp = cp, endogAA = endogAA)

end

dec = solveHHproblemEGM(param, 0.04, 1.0)

function solveInvariant(param, decisions)
	@unpack nS, nA, gA, transS = param
	@unpack ap, = decisions

	# Define aux parameters:
	maxiterInv = 50000
	tolInv = 1E-10

	# === 1. RETRIEVE GRID AND WEIGHT OF INTERPOLATED POLICIES ============== #
	# ps. the interpolation assume that policy is monotone
	ibelow = fill(0, nA, nS)
	iweight = zeros(nA, nS)
	for is = 1:nS
		iweight[:, is], ibelow[:, is] = coordGridIntp!(gA, ap[:, is], ibelow[:, is], iweight[:, is], robust = true)
	end
	# iweight is probability agent ends in grid "ibelow". 

	# ================ 2. ITERATE FORWARD DISTRIBUTION ======================== #
	dsn = fill(1.0/(nS*nA), (nA, nS)) # initial guess, positive mass everywhere must sum to one

	# Iterate over invariate distribution
    for iter = 1:maxiterInv

		# compute next distribution
		dsnNew = zeros(nA, nS)
		for ia = 1:nA, is = 1:nS
			if dsn[ia, is] > 0.0 # this speed a bit in some problems
				dsnNew[ibelow[ia, is], is] += iweight[ia, is] * dsn[ia, is]
				dsnNew[ibelow[ia, is] + 1, is] += (1-iweight[ia, is]) * dsn[ia, is]
			end
		end
		dsnNew = dsnNew*transS # apply the markov-chain of labor process

		# ===== Check if distribution has converged =====
        d = maximum(abs.(dsn - dsnNew))
		dsn = dsnNew # update distribution

        if d < tolInv
			println("Tol. achieved: $d")
			break
		end

		if iter == maxiterInv;
        	println("Max iterations achieved. Invariant distribution did not converge")
        end
		#println("Iter: $iter")
		#println("Tol: $d") # ps these prints will not show in the notebook.
	end

	return dsn

end

# Almost no one with the max asset grid. Around 1% of agents are constrained!
begin
    dsn1 = solveInvariant(param, dec)
    dsn2 = fill(1.0/(7*300), (300, 7)) 
        
        A = sum(dsn1.*repeat(param.gA, 1, param.nS)) # Aggregate asset supply
    #repeat(param.gA, 1, param.nS)
        #sum(dsn1[1,:]) # a_1 = 0.0
        #sum(dsn1)
end

## ======================================================================== ##
# This function solves one iteration of the model given prices
## ======================================================================== ##

function ExcessDemand(param, r)

	# Compute Capital Demand and implied wage from the production function
	@unpack Lbar, gA, nS, α, δ = param

	Kd = (α/(r+δ))^(1/(1-α))*Lbar # capital demand
	w = (1-α)*(Kd/Lbar)^α

	# Solve HH problem
	decisions = solveHHproblemEGM(param, r, w)

	# Compute Invariant distribution
	dsn = solveInvariant(param, decisions)

	# Compute Excess Demand of asset market
	Ea = sum(dsn.*repeat(gA, 1, nS)) # Aggregate asset supply
	excDem = (Ea -Kd)/((Ea+Kd)/2) # excess demand in percentage 
	#excDem = (Ea -Kd)


    return (excDem, decisions, dsn, w, Kd, Ea)
end


function ModelSolution(param)

	# Parameters
	r0 = 0.001 #  lower bound guess (make sure excess demand is negative)
	r1 = 1/param.β - 1 # upper bound guess (make sure excess demand is positive)

	tolEq = 0.001

	function objFct(rguess);
		println("\nInterest Rate Guess: $rguess")
		(excDem, ) = ExcessDemand(param, rguess);
		println("\nExcess Demand: $excDem")
		return excDem
	end

	r = fzero(objFct, (r0, r1), Roots.Brent(), atol = tolEq)
	(excDem, decisions, dsn, w, Kd, Ea) =  ExcessDemand(param, r) # get the stuff from the model

	return (decisions, dsn, w, r, Kd, Ea) # model stats
end


(decisions, dsn, w, r, Kd, Ea) = ModelSolution(setPar())
begin
	# Some useful variables
	@unpack α, Lbar, nS, nA, gA, gS = param
	Y = Kd^α*Lbar^(1-α) # agg output
	asset_distribution = sum(dsn, dims = 2)[:] # sum in the labor dimension
	const_hh = asset_distribution[1] 		# share of constrained households:
	
	# function that calculates gini (see wikipedia formula):
	function calGini2(dist, values)
		cumulative = cumsum(dist.*values)/sum(dist.*values) # cum distribution in %
		B = sum(cumulative.*dist) # total area below lorenz curve
	
		# with cns distribution this should be 0.5, but since it is a histogram it may deviate a bit
		AreaBelow45 = sum(dist.*(cumsum(dist)./sum(dist))) # A + B
		gini = (AreaBelow45 - B)/AreaBelow45
		return gini
	end

	gini_wealth = calGini2(asset_distribution,gA)
	
	# income/consumption, we must sort all the values first
	dnsvec = dsn[:]
		
    gYY = w.*repeat(gS', nA, 1) .+ (1.0+r).*repeat(gA, 1, nS) # income matrix
	inc_not_sorted = gYY[:] # income values vector
	idx_inc = sortperm(inc_not_sorted) # index of income	
	gini_income = calGini2(dnsvec[idx_inc],inc_not_sorted[idx_inc])

	cp_not_sorted = decisions.cp[:] # consumption values vector
	idx_cp = sortperm(cp_not_sorted) # index of consumption
	gini_consump = calGini2(dnsvec[idx_cp],cp_not_sorted[idx_cp])
	

	with_terminal() do # this is just to print stuff in the notebook
		println("Eq. Interest Rate and Wages: $r $w")
		println("Aggregate Capital and Asset Supply: $Kd $Ea")
		println("Aggregate Output: $Y")
		println("Fraction of constrained households: $const_hh")
		println("Gini of Wealth: $gini_wealth")
		println("Gini of Income: $gini_income")
		println("Gini of Consumption: $gini_consump")

	end
end


begin
	plot(gA, [decisions.ap[:,1] decisions.ap[:,nS]], label=["S min" "S max"], lw = 1, legend=:bottomright)
	title!("Policy Function")
	axis_lim = 30
	xlims!(0,axis_lim) 
	ylims!(0,axis_lim) 
	xlabel!("a")
	ylabel!("g_a")
end


begin
	plot(param.gA, asset_distribution, label="Aggregate", lw = 3)
	title!("Invariant Distribution")
	xlims!(-1,50) # little mass over 50
	ylims!(0,0.05) # Ps the mass of constrained HH is about 20% of individuals in S)1
	xlabel!("a")
	ylabel!("density")
end


begin
	plot(param.gA, [dsn[:,1]./sum(dsn[:,1]) dsn[:,nS]./sum(dsn[:,nS])], label=["S min" "S max" ], lw = 3)
	title!("Invariant Distribution")
	xlims!(-1,50) # little mass over 50
	ylims!(0,0.05) # Ps the mass of constrained HH is about 20% of individuals in S)1
	xlabel!("a")
	ylabel!("density")
end