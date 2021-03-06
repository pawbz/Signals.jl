module DSP

using Interpolation
using Conv
using Distributions
using DistributedArrays
using DSP # from julia
using SparseArrays
using LinearAlgebra
using FFTW



function chirp!(x; tgrid=nothing,
	fmin=nothing, # minimum fraction of Nyquist sampling
	fmax=nothing, # maximum fraction of Nyquist
	)

	fs=inv(step(tgrid))
	freqmin=(fmin===nothing) ? 0.0 : fmin*fs*0.5
	freqmax=(fmax===nothing) ? fs*0.5 : fmax*fs*0.5

	k=(freqmax-freqmin)*inv(tgrid[end]-tgrid[1])

	for it in 1:length(tgrid)
		t=(it-1)*step(tgrid)
		x[it] = sin(2. * π *(freqmin*t+k/2.0*t*t))
	end

end
"""
Random source signal models for Param.
Generate a band-limited 
random signal of a particular maximum time, `tmax`.
Sinc function in the frequency domain corresponds to 
a box car function is applied in the time domain.

This method has some extra allocations.
"""
function freq_rand!(x;
	tgrid=nothing,
	fmin=nothing, # minimum fraction of Nyquist sampling
	fmax=nothing, # maximum fraction of Nyquist
	nevent_fracs::Array{Float64}=[0.0],
	rng=Normal(),
	verbose=false,
	)

	nt=size(x,1)
	(tgrid===nothing) && (tgrid=range(0.0, stop=(nt-1)*1.0, step=1.0))  # create dummy tgrid
	fgrid=DSP.rfftfreq(length(tgrid),inv(step(tgrid))) # corresponding fgrid
	(nt≠length(tgrid)) && error("dimension")

	tmax=abs(tgrid[nt]-tgrid[1])

	xfreq=complex.(zeros(length(fgrid))) # allocation in freq domain
	xx=zeros(length(tgrid)) # allocation in time domain

	fs=inv(step(tgrid))
	freqmin=(fmin===nothing) ? 0.0 : fmin*fs*0.5
	freqmax=(fmax===nothing) ? fs*0.5 : fmax*fs*0.5
	(freqmax ≤ freqmin) && error("fmin and fmax")

	Δf = inv(tmax) # resolution in the frequency domain
	(Δf ≤ step(fgrid)) ? error("desired resolution smaller than grid sampling") :
	(Δf ≥ (freqmax-freqmin)) ? error("need to increase tmax") :
	fvec = [f for f in freqmin:Δf:freqmax]
	ifvec = fill(0, size(fvec))

	verbose && println("number of frequencies added to signal:\t", size(fvec,1))
	verbose && println("interval between random variable:\t", Δf)
	verbose && println("minimum frequency added to signal\t",minimum(fvec))
	verbose && println("maximum frequency added to signal\t",maximum(fvec))
	for iff in eachindex(fvec)
		ifvec[iff] =  Interpolation.indminn(fgrid, fvec[iff])[1]
	end

	for iff in eachindex(fvec)
		# translated sinc function
		fshift=fgrid[ifvec[iff]]  # shift for translated sinc function
		X=rand(rng) * exp(im*rand(Uniform(-pi, pi))) # uniformly distributed phase and random amplitude
		Xe=complex(1.0) # intialization for events
		for tfrac in nevent_fracs
			Xe += exp(im*fshift*2*π*tfrac*tgrid[nt]) # phase translations for each event
		end
		X *= Xe # add events to X
		for ifff in 1:length(fgrid)
			α=tmax*(abs(fgrid[ifff] - fshift)) # argument of sinc function
			xfreq[ifff] += (X * sinc(α)) # go
		end
	end
	xx=irfft(xfreq, nt)
	normalize!(xx)
	copyto!(x, xx)
	return x

end


"""
Tapering is necessary to be able to input random signal into finite-difference code
Filtering tapering are applied only if the length of the time series is greater than 10
"""
function get_tapered_random_tmax_signal(tgrid; 
					fmin=nothing,
					fmax=nothing,
					tmaxfrac::Float64=1.0,
					dist=Uniform(-2.0, 2.0),
					sparsep=1.0,
					taperperc=20.
					)
	filt_flag=(length(tgrid) > 5) && (!(fmin === nothing)) && (!(fmax===nothing))

	fs = inv(step(tgrid));
	if(filt_flag)
		designmethod = Butterworth(6);
		filtsource = Bandpass(fmin, fmax; fs=fs);
	end

	itind = argmin(abs.(tgrid.-abs(tmaxfrac)*tgrid[end]))
	if(tmaxfrac>0.0)
		its=1:itind
	elseif(tmaxfrac<0.0)
		its=itind:length(tgrid)
	end
	# 20% taper window
	twin = taper(ones(length(its)),taperperc) 
	X = zeros(length(its))
	wavsrc = zeros(length(tgrid)) 
	if(filt_flag) 
		X[:] = rand(dist, length(its)) .* twin
	else
		X[:] = rand(dist, length(its))
	end
	if(sparsep ≠ 1.0)
		Xs=sprandn(length(X), sparsep)
		X[findall(Xs.==0.0)].=0.0
	end
	# band limit
	(filt_flag) && (filt!(X, digitalfilter(filtsource, designmethod), X))
	
	(length(X) ≠ 1) && normalize!(X)
	wavsrc[its] = X
	return wavsrc
end


"""
* `x` : first dimension is time
* `p` : sparsity
* `rng` : Random Number Generator
* `btperc` : Taper Perc
* `etperc` : Taper Perc
* `bfrac` : zeros at the beginning of each coloumn
* `efrac` : zeros fraction at the end of each coloumn
"""
function tapered_rand!(x;p=1.0,rng=Normal(),btperc=0., etperc=0.,bfrac=0.0, efrac=0.0)
	x[:]=0.0
	nt=size(x,1)

	a=1+round(Int,bfrac*nt)
	b=nt-round(Int,efrac*nt)

	dd=size(x)[2:end]
	for i in CartesianRange(dd)
		xx=view(x,a:b,i)
		rand!(rng, xx)
		taper!(xx, bperc=btperc, eperc=etperc)
		nxx=size(xx,1)
		if(p ≠ 1.0)
			xs=sprandn(nxx, p)
			xx[findn(xs.==0.0)]=0.0
		end
	end
	scale!(x, inv(vecnorm(x)))
	return x
end

"""
randomly shift and add x to itself
using circshift; not memory efficient
"""
function cshift_add_rand!(x::AbstractVector; ntimes=1, tminfrac=0.0, tmaxfrac=1.0, rand_amp_flag=true)
	nt=size(x,1)
	a=1+round(Int,tminfrac*nt)
	b=round(Int,tmaxfrac*nt)

	for it in 1:ntimes
		if(rand_amp_flag)
			amp=randn()
		else
			amp=1.0
		end
		# a random shift
		its=rand(DiscreteUniform(a,b))
		xx = circshift(x,(its,)) # pre-allocate for performance?
		for i in eachindex(x)
			x[i] += xx[i] * amp
		end
	end
end


"""
Construct Toy Green's functions
Decaying peaks, control number of events, and their positions, depending on bfrac and efrac.
* `afrac` : control amplitude of peaks
"""
function toy_green!(x;nevents=1,afrac=[inv(2^(ie-1)) for ie in 1:nevents],bfrac=0.0,efrac=0.0)
	nt=size(x,1)
	nr=size(x,2)
	x[:]=0.0
	a=1+round(Int,bfrac*nt)
	b=nt-round(Int,efrac*nt)
	(b-a<nevents+2) && error("not enough samples")

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


function findfreq(
		  x::Array{Float64, ND},
		  tgrid;
		  attrib::Symbol=:peak,
		  threshold::Float64=-50.
		  ) where {ND}

	cx=rfft(x,[1]);
	fgrid=DSP.rfftfreq(length(tgrid),inv(step(tgrid))) # corresponding fgrid

	ax = (abs.(cx).^2); # power spectrum in dB

	if(maximum(ax) == 0.0)
		warn("x is zero"); return 0.0
	else 
		ax /= maximum(ax);
		ax = 10. .* log10.(ax)
	end

	if(attrib == :max)
		iii=findlast(ax .>= threshold)
	elseif(attrib == :min)
		iii=findfirst(ax .>= threshold)
	elseif(attrib == :peak)
		iii=argmax(ax)
	end
	ii=CartesianIndices(size(ax))[iii][1]
	return fgrid[ii]

end


"""
Cosine taper a N-dimensional array along its first dimension.

# Arguments
* `x::Array{Float64,N}` : 
* `perc::Float64` : taper percentage
"""
function taper(x::AbstractArray,perc::Float64)
	xout=copy(x);
	taper!(xout,perc)
	return xout
end

function taper!(x, perc=0.0; bperc=perc,eperc=perc)

	nt=size(x,1)
	nttb=min(round(Int,nt*bperc*0.01), nt)
	ntte=min(round(Int,nt*eperc*0.01), nt)
	kb=inv(2.0*round(Int,nt*bperc*0.01)-1)*pi
	ke=inv(2.0*round(Int,nt*eperc*0.01)-1)*pi
	for i in CartesianIndices(size(x))
		if(1≤i[1]≤nttb)
			x[i] *= sin((i[1]-1)*kb)
		end
		if(nt-ntte+1≤i[1]≤nt)
			x[i] *= sin((-i[1]+nt)*ke)
		end
	end
end

include("Filters.jl")


end # module
