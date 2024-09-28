function coordGridIntp!(xgrid::Array{Float64,1}, xtarget::Array{Float64,1}, 		
    ibelow::Array{Int64,1}, iweight::Array{Float64,1}; robust = false)
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
