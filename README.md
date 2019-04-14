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
  
Ultimately, when it comes down to it, the allan variance / allan deviation is nothing more than a Bode plot that analyzes the magnitude of different noise sources with differing frequency power.  
  
Allan variance power-law response
Power-law noise type	Phase noise slope	Frequency noise slope	Power coefficient	Phase noise
{\displaystyle S_{x}(f)} {\displaystyle S_{x}(f)}	Allan variance
{\displaystyle \sigma _{y}^{2}(\tau )} {\displaystyle \sigma _{y}^{2}(\tau )}	Allan deviation
{\displaystyle \sigma _{y}(\tau )} {\displaystyle \sigma _{y}(\tau )}
white phase modulation (WPM)	{\displaystyle f^{0}=1} f^{0}=1	{\displaystyle f^{2}} f^{2}	{\displaystyle h_{2}} h_{2}	{\displaystyle {\frac {1}{(2\pi )^{2}}}h_{2}} {\displaystyle {\frac {1}{(2\pi )^{2}}}h_{2}}	{\displaystyle {\frac {3f_{H}}{4\pi ^{2}\tau ^{2}}}h_{2}} {\displaystyle {\frac {3f_{H}}{4\pi ^{2}\tau ^{2}}}h_{2}}	{\displaystyle {\frac {\sqrt {3f_{H}}}{2\pi \tau }}{\sqrt {h_{2}}}} {\displaystyle {\frac {\sqrt {3f_{H}}}{2\pi \tau }}{\sqrt {h_{2}}}}
flicker phase modulation (FPM)	{\displaystyle f^{-1}} f^{-1}	{\displaystyle f^{1}=f} f^{1}=f	{\displaystyle h_{1}} h_{1}	{\displaystyle {\frac {1}{(2\pi )^{2}f}}h_{1}} {\displaystyle {\frac {1}{(2\pi )^{2}f}}h_{1}}	{\displaystyle {\frac {3[\gamma +\ln(2\pi f_{H}\tau )]-\ln 2}{4\pi ^{2}\tau ^{2}}}h_{1}} {\displaystyle {\frac {3[\gamma +\ln(2\pi f_{H}\tau )]-\ln 2}{4\pi ^{2}\tau ^{2}}}h_{1}}	{\displaystyle {\frac {\sqrt {3[\gamma +\ln(2\pi f_{H}\tau )]-\ln 2}}{2\pi \tau }}{\sqrt {h_{1}}}} {\displaystyle {\frac {\sqrt {3[\gamma +\ln(2\pi f_{H}\tau )]-\ln 2}}{2\pi \tau }}{\sqrt {h_{1}}}}
white frequency modulation (WFM)	{\displaystyle f^{-2}} f^{{-2}}	{\displaystyle f^{0}=1} f^{0}=1	{\displaystyle h_{0}} h_{0}	{\displaystyle {\frac {1}{(2\pi )^{2}f^{2}}}h_{0}} {\displaystyle {\frac {1}{(2\pi )^{2}f^{2}}}h_{0}}	{\displaystyle {\frac {1}{2\tau }}h_{0}} {\displaystyle {\frac {1}{2\tau }}h_{0}}	{\displaystyle {\frac {1}{\sqrt {2\tau }}}{\sqrt {h_{0}}}} {\displaystyle {\frac {1}{\sqrt {2\tau }}}{\sqrt {h_{0}}}}
flicker frequency modulation (FFM)	{\displaystyle f^{-3}} f^{{-3}}	{\displaystyle f^{-1}} f^{-1}	{\displaystyle h_{-1}} h_{{-1}}	{\displaystyle {\frac {1}{(2\pi )^{2}f^{3}}}h_{-1}} {\displaystyle {\frac {1}{(2\pi )^{2}f^{3}}}h_{-1}}	{\displaystyle 2\ln(2)h_{-1}} {\displaystyle 2\ln(2)h_{-1}}	{\displaystyle {\sqrt {2\ln(2)}}{\sqrt {h_{-1}}}} {\displaystyle {\sqrt {2\ln(2)}}{\sqrt {h_{-1}}}}
random walk frequency modulation (RWFM)	{\displaystyle f^{-4}} f^{{-4}}	{\displaystyle f^{-2}} f^{{-2}}	{\displaystyle h_{-2}} h_{{-2}}	{\displaystyle {\frac {1}{(2\pi )^{2}f^{4}}}h_{-2}} {\displaystyle {\frac {1}{(2\pi )^{2}f^{4}}}h_{-2}}	{\displaystyle {\frac {2\pi ^{2}\tau }{3}}h_{-2}} {\displaystyle {\frac {2\pi ^{2}\tau }{3}}h_{-2}}	{\displaystyle {\frac {\pi {\sqrt {2\tau }}}{\sqrt {3}}}{\sqrt {h_{-2}}}} {\displaystyle {\frac {\pi {\sqrt {2\tau }}}{\sqrt {3}}}{\sqrt {h_{-2}}}}
