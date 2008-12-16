% Module WAVEMODELS in WAFO Toolbox. 
% Version 2.1.1   07-Sep-2005 
%
% Readme     - New features, bug fixes, and changes in WAVEMODELS.
%
%Probability density functions (PDF)
% b04jpdf    - Brodtkorb (2004) joint (Scf,Hd) PDF from Japan Sea.
% b04jpdf2   - Brodtkorb (2004) joint (Scf,Hd) PDF from Japan Sea.
% b04pdf     - Brodtkorb (2004) joint (Scf,Hd) PDF of laboratory data.
% b04pdf2    - Brodtkorb (2004) joint (Scf,Hd) PDF of laboratory data.
% bmr00pdf   - Brodtkorb et.al (2000) joint (Scf,Hd) PDF from North Sea.
% bmr00pdf2  - Brodtkorb et.al (2000) joint (Scf,Hd) PDF from North Sea.
% pdfbray     - Beta Rayleigh PDF of wave heigth
% cav76pdf   - Cavanie et al. (1976) approximation of the density  (Tc,Ac)
% jhspdf     - Joint (Scf,Hd) PDF for linear waves with JONSWAP spectra.
% jhspdf2    - Joint (Scf,Hd) PDF for linear waves with a JONSWAP spectrum.
% jhsnlpdf   - Joint (Scf,Hd) PDF for nonlinear waves with a JONSWAP spectra.
% jhsnlpdf2  - Joint (Scf,Hd) PDF for non-linear waves with a JONSWAP spectra.
% jhvnlpdf   - Joint (Vcf,Hd) PDF for linear waves with a JONSWAP spectrum.
% jhvnlpdf2  - Joint (Vcf,Hd) PDF for non-linear waves with a JONSWAP spectrum.
% jhvpdf     - Joint (Vcf,Hd) PDF for linear waves with a JONSWAP spectrum.
% jhvpdf2    - Joint (Vcf,Hd) PDF for linear waves with a JONSWAP spectrum. 
% lh83pdf    - Longuet-Higgins (1983) approximation of the density (Tc,Ac) 
% pdfseastate - Long Term Joint Sea State PDF of Hm0 and wave period 
% mk87pdf    - Myrhaug and Kjeldsen (1987) joint (Scf,Hd) PDF. 
% mk87pdf2   - Myrhaug and Kjeldsen (1987) joint (Scf,Hd) PDF. 
% ochi98pdf  - Ochi's (1998) PDF of peaks and troughs of non-gaussian processes 
% ohhpdf     - Marginal wave height, Hd, PDF for Bimodal Ochi-Hubble spectra. 
% ohhspdf    - Joint (Scf,Hd) PDF for linear waves with Ochi-Hubble spectra. 
% ohhspdf2   - Joint (Scf,Hd) PDF for linear waves with Ochi-Hubble spectra. 
% ohhsspdf   - Joint (Scf,Hd) PDF for linear waves in space with Ochi-Hubble spectra. 
% ohhsspdf2  - Joint (Scf,Hd) PDF for linear waves in space with Ochi-Hubble spectra. 
% ohhvpdf    - Joint (Vcf,Hd) PDF for linear waves with Ochi-Hubble spectra. 
% tay81pdf   - Tayfun (1981) PDF of  breaking limited wave heights 
% tay90pdf   - Tayfun (1990) PDF of large wave heights 
% thpdf      - Marginal wave height, Hd, pdf for Torsethaugen spectra. 
% thsnlpdf2  - Joint (Scf,Hd) PDF for nonlinear waves with Torsethaugen spectra.
% thsnlpdf   - Joint (Scf,Hd) PDF for nonlinear waves with Torsethaugen spectra.
% thspdf     - Joint (Scf,Hd) PDF for linear waves with Torsethaugen spectra. 
% thspdf2    - Joint (Scf,Hd) PDF for linear waves with Torsethaugen spectra. 
% thsspdf    - Joint (Scf,Hd) PDF for linear waves in space with Torsethaugen spectra. 
% thsspdf2   - Joint (Scf,Hd) PDF for linear waves in space with Torsethaugen spectra. 
% thvpdf     - Joint (Vcf,Hd) PDF for linear waves with Torsethaugen spectra. 
% thvpdf2    - Joint (Vcf,Hd) PDF for linear waves with Torsethaugen spectra. 
% trraylpdf  - Transformed Rayleigh approximation for amplitudes 
%
%Cumulative distribution functions (cdf)
% b04cdf     - Brodtkorb (2004) joint (Scf,Hd) CDF of laboratory data.
% b04jcdf    - Brodtkorb (2004) joint (Scf,Hd) CDF from Japan Sea.
% bmr00cdf   - Brodtkorb et.al (2000) joint (Scf,Hd) CDF from North Sea.
% cdfbray     - Beta Rayleigh CDF of wave heights 
% jhscdf     - Joint (Scf,Hd) CDF for linear waves with a JONSWAP spectrum.
% jhsnlcdf   - Joint (Scf,Hd) CDF for non-linear waves with a JONSWAP spectrum.
% jhvcdf     - Joint (Vcf,Hd) CDF for linear waves with a JONSWAP spectrum. 
% jhvnlcdf   - Joint (Vcf,Hd) CDF for non-linear waves with a JONSWAP spectrum. 
% mk87cdf    - Myrhaug and Kjeldsen (1987) joint (Scf,Hd) CDF. 
% ochi98cdf  - Ochi's (1998) CDF of peaks and troughs of non-gaussian processes 
% ohhcdf     - Marginal wave height, Hd, CDF for Ochi-Hubble spectra.
% ohhscdf    - Joint (Scf,Hd) CDF for linear waves with Ochi-Hubble spectra.
% ohhsscdf   - Joint (Scf,Hd) CDF for linear waves in space with Ochi-Hubble spectra.
% tay81cdf   - Tayfun (1981) CDF of  breaking limited wave heights 
% tay90cdf   - Tayfun (1990) CDF of large wave heights 
% thcdf      - Marginal wave height, Hd, cdf for Torsethaugen spectra. 
% thscdf     - Joint (Scf,Hd) CDF for linear waves with Torsethaugen spectra. 
% thsnlcdf   - Joint (Scf,Hd) CDF for non-linear waves with Torsethaugen spectra.
% thsscdf    - Joint (Scf,Hd) CDF for linear waves in space with Torsethaugen spectra.
% thvcdf     - Joint (Vcf,Hd) CDF for linear waves with Torsethaugen spectra. 
%
%Random number generators
% mk87rnd    - Random points from MK87 distribution of steepness and wave height.
%
%Parameter estimation 
% fitbray     - Parameter estimates for Beta-Rayleigh data. 
% ochi98fit  - Parameter estimates and confidence intervals for Ochi data. 
%
%Misc
% jhwparfun  - Wave height, Hd, distribution parameters for Jonswap spectrum. 
% jhnlwparfun - Wave height, Hd, distribution parameters for Stokes waves with Jonswap spectrum.
% ohhgparfun - Wave height, Hd, distribution parameters for Ochi-Hubble spectra. 
% tay81fun   - Internal function to tay81pdf and tay81cdf. 
% tay90fun   - Internal integrand function to tay90pdf 
% thgparfun  - Wave height, Hd, distribution parameters for Torsethaugen spectra. 
% thwparfun  - Wave height, Hd, distribution parameters for Torsethaugen spectra.
%
%
 



 





