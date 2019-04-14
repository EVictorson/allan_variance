# allan_variance
Synthetic IMU error models based on first order Gauss-Markov process and corresponding allan variance analysis
  
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
