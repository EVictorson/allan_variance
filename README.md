# allan_variance
Synthetic IMU error models based on first order Gauss-Markov processes and corresponding allan variance analysis.
  
I'm hoping to expand upon this at some point in the future and create more of an educational tool, but for now here is just a dump of things I have learned.

educational sources:  
https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=660628  IEEE std 953-1997  
http://home.engineering.iastate.edu/~shermanp/AERE432/lectures/Rate%20Gyros/14-xvagne04.pdf  
http://cache.freescale.com/files/sensors/doc/app_note/AN5087.pdf  
https://openi.nlm.nih.gov/detailedresult.php?img=PMC3812568_sensors-13-09549f4&req=4  
https://etd.ohiolink.edu/!etd.send_file?accession=osu1420709961&disposition=inline#page=144&zoom=100,0,661  
https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=494457 IEEE std 647-1995 (bad qc and tc calculation)  
note: when I say bad qc and tc calculation, I mean I believe IEEE to have published a mistake in this  
version of the standard.  
  
Ultimately, when it comes down to it, the allan variance / allan deviation is nothing more than a Bode plot that analyzes the magnitude of different noise sources with differing frequency power.  This will cause different contributing noise terms to appear with differing slopes on the allan deviation plot.  
  
Quantization Noise: slope -1  
Angle Random Walk: slope = -1/2  
Bias Instability: slope = 0  
Rate Random Walk: slope = 1/2  
Rate Ramp: slope = 1  
  
From the navigation class part 2:
Nav_class part 2 section 5-9
As previously stated, the gyro bias error is typically modeled as a "slow" 
Gauss-Markov process (1000 hour correlation time) which can be considered 
a random constant for a single flight.  The typical gyro bias error budget
is on the order of 0.01 deg/hr.
   
Nav class part 2 section 5-11
As stated before, the gyro scale factor error is also typically modeled as
a Gauss-Markov process with a long correlation time (1000 hours).  Again, 
for a single flight this can be considered as a random constant. 
The typical gyro scale factor error budget is on the order of 5 
parts-per-million RMS
  
Again, as with the gyro, the accel bias error is typically modeled as a 
"slow" Gauss-Markov process which can be considered a random constant for 
a single flight.  The typical accel bias error budget is on the order 
of 84 µg.
   
The accel scale factor error is also typically modeled as a Gauss-Markov 
process with a long correlation time (100 hours).  Again, for a single 
flight this can be considered as a random constant. 
The typical accel scale factor error budget is on the order of 300 
parts-per-million RMS.
  
The following terms are neglected, but documented for future use if
necessary: scale factor nonlinearity, nonorthogonality, misalignment
  
Validation of this class has been performed via allan deviation
analysis, wherein each individual process (bias, bias inrun, scale factor)
was analyzed and the time constant (Tc) and noise magnitude (qc) were
verified.  For performing this type of analysis again see the
allan_variance_testing.m file.
  
At its core, an exponentially correlated, noise driven stochastic
process (gauss markov process) may be modeled with the following
differential equation: dx/dt = (-1/tau) * x + N(mu, sigma)
which is just an exponentially decaying autoregressive term with
driven white noise.
  
The discrete time version of this is:   
  
x(n+1) = x(n) * exp((-1/tau) * dt) + N(mu,sigma) * dt
   
The first order gauss-markov process may also be defined in terms of
its autocorrelation function:   
  
Rxx = (sigma^2) * exp(-beta * tau)  
  
where sigma is the ENTIRE process variation (not the driven noise variation), beta is
the inverse of the correlation time T = 1/beta (time constant), and tau is
the autocorrelation time lag.  For sufficiently long time series (much 
longer than the time constant), the autocorrelation plot of the first order 
Gauss-Markov process can verify the correlation time, which is the time 
lag at which Rxx = (sigma^2) / e, and the entire time series variance 
will appear as the peak at time lag = 0.
Note that the driven noise standard deviation cannot be validated via
autocorrelation. Also note that matlab has two different autocorrelation
functions that will produce slightly different results if your time scales
are short (autocorr and xcorr), as well as an autocovariance function, xcov.
Autocorr and xcov are normalized by the process mean such that a constant
process will have a steady autocorrelation (autocorrelation is just
autocovariance scaled by the inverse of the process variance). xcorr,
however, does not do this, so a constant process will have a decaying
autocorrelation as the time lag increases, which results in a shortening
of the summation by 1 element each iteration.
     
For process lengths on the order of magnitude of the time constant or
less the process will be dominated by the integrated white noise term, resulting
in a random walk.  For time scales much longer than the time constant the
exponentially correlated term will begin to be apparent.  This can be see
on an allan deviation plot, where for sampling intervals much shorter than
the time constant the gauss-markov allan variance reduces to that of a
singly integrated white noise process (rate random walk), whose slope 
is +1/2, and the noise magnitude (standard deviation) may be picked off 
by finding the intersection of the +1/2 slope line and sampling_interal = 3.  
This can also be verified by differentiating this time series to obtain 
gaussian noise, whose allan deviation has a slope of -1/2, and the noise
magnitude is interpreted as the intersection of thise line with
sampling_interval = 1.
  
For process time scales large enough to see the decay time constant, and
allan deviations with sampling intervals ranging from much less than the
time constant to sampling intervals much larger than the time constant,
the slope will transition from +1/2 for sampling intervals much shorter 
than the time constant, to 0 as the sampling interval approaches
the time constant at it's maxima, to -1/2 as the interval becomes much
larger than the time constant.
   
When multiple gauss markov processes are added together it will be
nearly impossible to pick apart any of these parameters, but when run
individually the driven noise magnitude (qc) and time constant (Tc) can be
found on the allan deviation plot as follows:  
Tc = argmax_sampling_interval / 1.89, where argmax_sampling_interval is
the allan deviation sampling interval that maximizes the allan deviation.
The driven noise magnitude may then be found with the following
relationship:  
  
qc = sigma_max / (0.437 * sqrt(Tc)) 
  
where sigma_max is the allan
deviation maxima.
  
If the combined output of all of these Gauss-Markov processes are
analyzed, the only easy thing to pick off will be the velocity random
walk / angle random walk, which may be found by fitting a line to the
segment with a slope of -1/2 and finding the intersection of this line
with tau = 1.
