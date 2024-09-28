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
include("coordGridIntp.jl")
param = setPar() # check the parameters above
include("solver.jl")

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