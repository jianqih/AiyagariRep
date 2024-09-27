function setPar(;	
    nA = 200, # Asset grid size
	nS = 2,		# Labor endowment grid size
	α = 0.33,
	δ = 0.1,
	β = 0.96,
	γ = 2.0,	# CRRA parameter
    ϕ = 0.0,   # borrowing constraint
	ρ = 0.9, # autocorrelation
	σ = 0.1, # std deviation

	amax = 30.0, # maximum grid point
	grid_growth = 0.00  # asset grid curvature (growth rate between points)

	)

	#### Define the labor grid
	# gS, transS = tauchen(nS, ρ, σ)
    # gS = exp.(gS)
    gS = [0.8; 1.2]
    transS =[0.9 0.1;0.1 0.9]

	# compute invariant distribution of labor process (by iteration)
	invS = ones(nS) / nS # initial guess
	tol=1E-11; maxit=10^4
	for it in 1:maxit
    	invS_new = (invS' * transS)'
    	if maximum(abs.(invS_new .- invS)) < tol; break; end
    	invS .= invS_new
	end

	#invS = transS^200 # compute invariant labor distribution
	Lbar = sum(gS.*invS) # aggregate labor supply

	#### Define the asset grid:
	if grid_growth == 0.0
		gA = collect(range(-ϕ, amax, length = nA)) # evenly sspaced grid
	elseif grid_growth>0.0
		gA = fill(0.0, nA) # pre-allocate grid
    	for i in 1:nA
    		gA[i] = -ϕ + (amax - (-ϕ) )*((1 + grid_growth)^(i-1.0) -1)/((1 + grid_growth)^(nA-1.0) -1)
    	end

	end

	return (α = α, δ = δ, β = β, γ = γ, ϕ = ϕ,
	nA = nA, gA = gA, nS = nS, gS = gS, transS = transS, Lbar = Lbar, invS = invS )
end
