%Module SPEC in WAFO Toolbox. 
%Version 2.5.2   07-Feb-2011 
% 
%   Readme               - Readme file for module SPEC in WAFO Toolbox 
%
%Model spectra
%   bretschneider        - Calculates (and plots) a Bretschneider spectral density.
%   demospec             - Loads a precreated spectrum of chosen type
%   getjonswappeakedness - Peakedness factor Gamma given Hm0 and Tp for JONSWAP
%   getjonswapseastate   - Return estimated seastate from windspeed and fetch
%   jonswap              - Calculates (and plots) a JONSWAP spectral density
%   mccormick            - Calculates (and plots) a McCormick spectral density.
%   mkbretschneider      - Function handle to BRETSCHNEIDER spectral density
%   mkjonswap            - Function handle to JONSWAP spectral density
%   mkochihubble         - Function handle to bimodal spectral density.
%   mkspreading          - Function handle to Directional spreading function.
%   mktmaspec            - Function handle to JONSWAP spectral density modified for finite water depth
%   mktorsethaugen       - Function handle to double peaked (swell + wind) spectrum 
%   ohspec               - Calculates (and plots) a Ochi-Hubble spectral density.
%   ohspec2              - Calculates and plots a bimodal Ochi-Hubble spectral density
%   ohspec3              - Calculates Bimodal Ochi-Hubble spectral densities
%   oscspec              - Spectral density for a harmonic oscillator
%   phi1                 - factor for transforming spectra to finite water depth spectra
%   pmspec               - Calculates (and plots) a Pierson-Moskowitz spectral density.
%   tmaspec              - Calculates a JONSWAP spectral density for finite water depth
%   torsethaugen         - Calculates a double peaked (swell + wind) spectrum 
%   wallop               - Calculates (and plots) a Wallop spectral density.
%
%Directional spectra and spreading functions
%   mkdspec              - Make a directional spectrum
%   spreading            - Directional spreading functions
%
%Spectral characteristics
%   dspec2char           - Evaluates directional spectral characteristics 
%   spec2bw              - Evaluates some spectral bandwidth and irregularity factors
%   spec2char            - Evaluates spectral characteristics and their covariance
%   spec2mom             - Calculates spectral moments from spectrum
%
%Spectrum and covariance functions
%   cov2spec             - Computes spectral density given the auto covariance function 
%   covinterp            - Interpolation of covariance function and derivatives
%   createcov            - Covariance class constructor
%   createspec           - Spectrum structure constructor
%   freqtype             - returns the frequency type of a Spectral density struct.
%   lagtype              - Returns the lag type of a Covariance struct.
%   rotspec              - Rotate spectrum clockwise around the origin. 
%   scalespec            - Scale spectral density so that the moments equals m0,m2.  
%   spec2cov             - Computes covariance function and its derivatives 
%   spec2cov2            - Computes covariance function and its derivatives, alternative version
%   spec2dt              - Computes sampling interval from Nyquist frequency in spectrum
%   spec2spec            - Transforms between different types of spectra
%   specinterp           - Interpolation and zero-padding of spectrum
%   ttspec               - Toggle Transform between angular frequency and frequency spectrum
%   normspec             - Normalize a spectral density such that m0=m2=1
%   specplot             - Plot a spectral density 
%
%Dispersion relation
%   k2w                  - Translates from wave number to frequency
%   w2k                  - Translates from frequency to wave number






