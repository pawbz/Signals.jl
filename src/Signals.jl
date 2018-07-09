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




"""
Construct Toy Green's functions
Decaying peaks, control number of events, and their positions, depending on bfrac and efrac.
* `afrac` : control amplitude of peaks
"""
function toy_reflec_green!(x; bfrac=0.0, c=1.0, afrac=1.0)
	nt=size(x,1)
	nr=size(x,2)
	c=c*sqrt((nr^2)/(nt-1))

	it0=1+round(Int,bfrac*nt)
	for ir in 1:nr
		it=round(Int, it0+(ir*ir)*inv(c)*inv(c))
		if(it<=nt)
			x[it, ir]=1.0*afrac
		end
	end
	return x
end

function toy_direct_green!(x; bfrac=0.0, c=1.0, afrac=1.0)
	nt=size(x,1)
	nr=size(x,2)
	c=c*((nr)/(nt-1))

	it0=1+round(Int,bfrac*nt)
	for ir in 1:nr
		it=round(Int, it0+(ir)*inv(c))
		if(it<=nt)
			x[it, ir]=1.0*afrac
		end
	end
	return x

end

#=

	a=1+round(Int,bfrac*nt)
	b=nt-round(Int,efrac*nt)
	(b-a<nevents+2) && error("not enough samples")

	mult_shift=zeros(Int64,nevents+1)
	itvec=zeros(Int64,nevents+1)
	for ir in 1:nr
		itvec[:]=0
		for ie in 1:nevents
			if(ie==1)
				# first event, direct arrival
				itvec[ie]=rand(DiscreteUniform(a,b-nevents-1))
			else
				# all other events after first event
				itvec[ie]=rand(DiscreteUniform(itvec[1]+1,b-1))
				if(mult_shift[ie]==0)
					mult_shift[ie]=rand(DiscreteUniform(1,b-1-itvec[1]))
				end

			end
			sgn=1.
			# randomly choose sign for all other events
			if(ie≠1)
				sgn=(rand(Bernoulli())==1) ? 1. : -1.
			end
			x[itvec[ie],ir]=sgn*afrac[ie]
		end
	end
	return x
end

=#

end # module
