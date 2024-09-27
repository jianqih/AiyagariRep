"""
    linearGridIntp!(xgrid, ygrid, xtarget, ytarget)

Returns a linear interpolation of the vector xtarget on grid (xgrid, ygrid).
Assume xgrid and xtarget is sorted.
# Arguments
- `xgrid::Array{AbstractFloat,1}` : xgrid.
- `ygrid::Array{AbstractFloat,1}` : ygrid.
- `xtarget::Array{AbstractFloat,1}` : grid to interpolate.
# Output (inplace)
- `ytarget::Array{AbstractFloat,1}`
"""
function linearGridIntp!(xgrid::Array{Float64,1}, ygrid::Array{Float64,1}, xtarget::Array{Float64,1}, ytarget::Array{Float64,1})

	# check if is sorted
	@assert (issorted(xgrid) & issorted(xtarget))

	xg = length(xgrid); #xt = length(xtarget)
	xi  = 1; xlow = xgrid[1]; xhigh = xgrid[2] # initialize
	@inbounds for (it, xtval) in enumerate(xtarget)
		while xi < xg - 1 # find grid point
			if xhigh >= xtval; break; end
			xi += 1
			xlow = xhigh
            xhigh = xgrid[xi + 1]
		end
		xprob = (xhigh - xtval) / (xhigh - xlow)
		ytarget[it] = ygrid[xi] * xprob + (1.0 - xprob) * ygrid[xi+1]
	end
	return ytarget
end