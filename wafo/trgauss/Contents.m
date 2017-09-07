% Module TRGAUSS in WAFO Toolbox. 
% Version 2.5.3   30-May-2017
% 
%   Readme         - New features, bug fixes, and changes in TRGAUSS.
%
% Misc
%   createpdf      - PDF struct constructor
%   pdfplot        - Plot contents of pdf structures
%   trplot         - Plots transformation, g, eg. estimated with dat2tr.
%
%Transforms and non-linearities
%   dat2gaus       - Transforms  x  using the transformation  g.
%   gaus2dat       - Transforms  xx  using the inverse of transformation  g.
%   testgaussian   - Test if a stochastic process is Gaussian.
%   spec2skew      - Estimates the moments of 2'nd order non-linear waves 
%   trangood       - Makes a transformation that is suitable for efficient transforms.
%   tranproc       - Transforms process X and up to four derivatives 
%   trmak          - Put together a transformation object.
%   troptset       - Create or alter TRANSFORM OPTIONS structure.
%   trunmak        - Split a transformation object into its pieces.
%
% Transformed Gaussian model estimation
%   cdf2tr         - Estimate transformation, g, from observed CDF.
%   dat2tr         - Estimate transformation, g, from data.
%   hermitetr      - Estimate transformation, g, from the first 4 moments.
%   ochitr         - Estimate transformation, g, from the first 3 moments.
%   lc2tr          - Estimate transformation, g, from observed crossing intensity.
%   lc2tr2         - Estimate transformation, g, from observed crossing intensity, version2.
%
%Gaussian probabilities and expectations
%   cdfnorm2d      - Bivariate Normal cumulative distribution function
%   prbnorm2d      - Bivariate Normal probability   
%   prbnormnd      - Multivariate Normal probability by Genz' algorithm.
%   prbnormndpc    - Multivariate Normal probabilities with product correlation
%   prbnormtnd     - Multivariate normal or T probability by Genz' algorithm.
%   prbnormtndpc   - Multivariate normal or T probability with product correlation structure.
%   rind           - Computes multivariate normal expectations
%   rindoptset     - Create or alter RIND OPTIONS structure.
%
%Probability density functions (pdf) or intensity matrices
%   chitwo2lc_sorm - SORM-approximation of crossing intensity, noncentral Chi^2 process 
%   chitwo2lc_sp   - Saddlepoint approximation of crossing intensity, noncentral Chi^2 process 
%   dirsp2chitwo   - Parameters in non-central CHI-TWO process for directional Stokes waves. 
%   iter           - Calculates a Markov matrix given a rainflow matrix
%   iter_mc        - Calculates a kernel of a MC given a rainflow matrix
%   mc2rfc         - Calculates a rainflow matrix given a Markov chain with kernel f_xy;
%   mctp2rfc       - Rainflow matrix given a Markov matrix of a Markov chain of turning points
%   mctp2tc        - Calculates frequencies for the  upcrossing troughs and crests
%   nt2fr          - Calculates the frequency matrix given the counting distribution matrix.
%   spec2cmat      - Joint intensity matrix for cycles (max,min)-, rainflow- and (crest,trough)
%   spec2mmtpdf    - Joint density of Maximum, minimum and period.
%   spec2tccpdf    - Density of crest-to-crest wave -period or -length.  
%   spec2thpdf     - Joint density of amplitude and period/wave-length characteristics
%   spec2tpdf      - Density of crest/trough- period or length.      
%   spec2tpdf2     - Density of crest/trough- period or length, version 2.     
%   specq2lc       - Saddlepoint approximation of crossing intensity for quadratic sea.
%   th2vhpdf       - Transform joint T-H density to V-H density
%
%Cumulative distribution functions (cdf)
%   cdflomax       - CDF for local maxima for a zero-mean Gaussian process
%   spec2AcAt      - Survival function for crests and troughs, R(h1,h2)=P(Ac>h1,At>h2).
%   spec2Acdf      - CDF for crests P(Ac<=h) or troughs P(At<=h).




