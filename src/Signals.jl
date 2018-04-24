module Signals

include("DSP.jl")
include("Wavelets.jl")


"""

Output xN after adding noise to x

SNR:: in decibel scale
"""
function addnoise!(xN, x, SNR)


	# store zero mean version of x in xN
	xavg=0.0
	for i in eachindex(x)
		xavg += x[i]
	end
	xavg /= length(x)

	for i in eachindex(x)
		xN[i] = x[i] - xavg
	end

	σx=var(xN)

	σxN=sqrt(σx^2*inv(10^(SNR/10.)))
	
	# factor to be multiplied to each scalar
	α=sqrt(σxN)
	for i in eachindex(x)
		xN[i] = x[i] + α*randn()
	end
end




end # module
