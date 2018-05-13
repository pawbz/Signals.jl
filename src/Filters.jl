
mutable struct P_filt{T<:Real}
	fs::Float64
	freq::Vector{Float64}
	p_conv::Conv.P_conv{T,1,1,1}
end

function P_filt(nt::Int=128, filt_func=x->butterworth(x, order=2, fmin=0.0, fmax=1.0); fs=1)
	# create a plan for conv
	p_conv=Conv.P_conv(gsize=[nt], dsize=[nt], ssize=[nt],)

	freq=rfftfreq(p_conv.np2, fs)
	filt=filt_func(freq)
	for i in eachindex(p_conv.sfreq)
		p_conv.sfreq[i] = filt[i]
	end


	A_mul_B!(p_conv.spad, p_conv.sifftp, p_conv.sfreq)
	Conv.pad_truncate!(p_conv.s, p_conv.spad, p_conv.slags[1], 
		    p_conv.slags[2], p_conv.np2, -1)

	fs=Float64(fs)
	return P_filt(fs, Vector(freq), p_conv)

end

function filt!(xx, x, pa::P_filt)

	for ii in 1:size(xx,2)
		xv=view(x, :, ii)
		xxv=view(xx, :, ii)
		Conv.mod!(pa.p_conv, :d, d=xxv, g=xv)
	end
end

"""
A bandpass butterworth filter using `fmin` and `fmax` in the frequency domain.
Return either a zero-phase or minimum-phase filter using `attrib`.
More info on minimum-phase filtering: http://www.katjaas.nl/minimumphase/minimumphase.html.
! then it is converted to minimum phase filter in the cepstral domain
! more info: http://www.katjaas.nl/minimumphase/minimumphase.html
! positive quefrencies correspond to minimum phase component of the signal
! negetive quefrencies correspond to maximum phase component of the signal
! signal = minimum phase [convolved with] maximum phase


# Arguments

* ``

# Keyword Arguments

* `order::Int64` : order of the butterworth filter
* `attrib::Symbol`
  * `=:zp` means zero phase 
  * `=:mp` means zero phase 
* `fmin::Float64`
  * `=0.` means a lowpass filter at fmax
* `fmax::Float64`
  * `=0.` means a highpass filter at fmin

"""
function butterworth(fvec; order::Int64=2, attrib::Symbol=:zp, fmin=0.0, fmax=0.0)
# Author : Pawan Bharadwaj
#          p.b.pisupati@tudelft.nl
# August 2017, imported to Julia from FORTRAN90
	(fmin < 0.0) && error("fmin cannot be .lt. zero")
	(fmax < 0.0) && error("fmax cannot be .lt. zero")
	(fmax ≠ 0.0) && (fmax <= fmin) && error("fmax .le. fmin") 



	if(fmin == 0.0) 
		# lowpass filter
		F  = complex.((1.0 + (fvec./fmax).^(2.0*order)).^(-1.0), 0.0)
	elseif(fmax == 0.0) 
		# highpass filter
		F  = complex.((1.0 + (fmin./fvec).^(2.0*order)).^(-1.0), 0.0)
	elseif((fmin ≠ 0.0) & (fmax ≠ 0.0))
		# bandpass filter
		x1 = sqrt(fmin*fmax) / (fmax-fmin) * (fvec./sqrt(fmin*fmax) + sqrt(fmax*fmin)./fvec);
		F  = complex.((1.0 + (x1).^(2.0*order)).^(-1.0), 0.0)
	end
	# DC component is always forced zero
	F[1] = complex(0.0, 0.0);

	# conversion to minimum phase
	if(attrib == :mp)
	# to prevent log(0)
		damp = 1e-20 * maximum(abs.(F))
		# logarithm 
		X = log.(complex.(abs.(F) + damp, 0.0)) 
		# to cepstral domain - IFFT
		ifft!(X)
		# only real part output
		X = complex.(real.(X), 0.0);
		# scaling
		# X = X / complex(real(npow2), 0.0)

		# positive cepstrum x 2
		X[2 : npow2/2 + 1] *= complex(2.0, 0.0)
		# remove negative quefrencies
		X[npow2/2 + 2 : npow2] = complex(0.0, 0.0) 

		# FFT
		fft!(X)

		# exponential
		F = exp.(X)

		F /= complex(maximum(abs.(F)), 0.0)
	end

	return F

end # butterworth


