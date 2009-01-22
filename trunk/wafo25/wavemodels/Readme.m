%README New features, bug fixes, and changes in WAVEMODELS.
% Version 2.1.1   07-Nov-2004
%
%  FIXES
%  ~~~~~
%  ochi98fit -  replaced call to fmins with fminsearch  
%  Contents.m file missed reference to some functions
%  The distributions:
%    thsspdf, thspdf, thsnlpdf, jhspdf and jhsnlpdf
%  are fitted to new simulations saved in:
%     jhsnlpar21-Jul-2004.mat, jhspar13-Jul-2004.mat 
%     thsnlpar06-Aug-2004.mat, thspar19-Jul-2004.mat
%     thsspar27-Jul-2004.mat.
%
%  NEW FUNCTIONS
% ~~~~~~~~~~~~~~~~
% b04jpdf    - Brodtkorb (2004) joint (Scf,Hd) PDF from Japan Sea.
% b04jpdf2   - Brodtkorb (2004) joint (Scf,Hd) PDF from Japan Sea.
% b04pdf     - Brodtkorb (2004) joint (Scf,Hd) PDF of laboratory data.
% b04pdf2    - Brodtkorb (2004) joint (Scf,Hd) PDF of laboratory data.
% bmr00pdf   - Brodtkorb et.al (2000) joint (Scf,Hd) PDF from North Sea.
% bmr00pdf2  - Brodtkorb et.al (2000) joint (Scf,Hd) PDF from North Sea.
%
% b04cdf     - Brodtkorb (2004) joint (Scf,Hd) CDF of laboratory data.
% b04jcdf    - Brodtkorb (2004) joint (Scf,Hd) CDF from Japan Sea.
% bmr00cdf   - Brodtkorb et.al (2000) joint (Scf,Hd) CDF from North Sea.
% 
% 
%---------------------------------------------------
% Version 2.1.0   07-Apr-2004
%
%  FIXES
%  ~~~~~
%  trraylpdf - added example
%
%  For the following files the input is changed and the help header is
%  updated accordingly:
%
%  mk87cdf 
%  mk87pdf 
%  mk87pdf2
%  ohhsspdf
%  ohhspdf
%  thsspdf 
%  thspdf
%  thvdf
%
%  NEW FUNCTIONS
% ~~~~~~~~~~~~~~~~
%
%   jhvnlpdf
%   ohsspdf2
%   ohhspdf2 
%   thvpdf2
%   thsspdf2
%   thspdf2 
%   ohhsscdf
%   ohhscdf
%   thvcdf
%   thsscdf
%   thscdf
%   
%   
%-------------------------------------------------------------------
% Version 2.0.5   14-Jun-2003 
% 
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
% 
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO/wavemodels Toolbox.
%
% List of changes:
% 
% NOTE:
% Items marked with:
%       (N)  : needs a new  name 
%        *   : not yet available
%        +   : needs more work
%        **  : have changed in a way which might affect your code!
%
% 
%  FIXES
%  ~~~~~
%   - cav76pdf : 1) fixed transform and normalisation using Lindgren &
%                   Rychlik (1982) paper. 
%                2)  fixed xlabels i.e. f.labx={'Tc','Ac'}
%
%
%  ENHANCEMENTS
%  ~~~~~~~~~~~~
%
%  The dependence on the MATLAB STATS toolbox is now removed. Only
%  functions from the WAFO/wstats directory is needed. All the loglikelihood 
%  functions are removed, due to the fact that they are no longer needed by 
%  the XXXXXfit.m functions.
%
%  NEW FUNCTIONS
% ~~~~~~~~~~~~~~~~
%
%  OHHSPDF   - Joint (Scf,Hd) pdf in time for bimodal Ochi-Hubble spectra.
%  OHHPDF    - Marginal wave height, Hd, pdf for Bimodal Ochi-Hubble spectra.
%  OHHVPDF   - Joint (Vcf,Hd) pdf for bimodal Ochi-Hubble spectra.
%  OHHSSPDF  - Joint (Scf,Hd) pdf in space for bimodal Ochi-Hubble spectra.
%  OHHGPARFUN - Wave height, Hd, distribution parameters for Ochi-Hubble spectra.
%  OHHCDF    - Marginal wave height, Hd, cdf for Bimodal Ochi-Hubble spectra.
%  THSSPDF   - Joint (Scf,Hd) pdf in space for Torsethaugen spectra.
%  THSPDF    - Joint (Scf,Hd) pdf for Torsethaugen spectra.
%  THPDF     - Marginal wave height, Hd, pdf for Torsethaugen spectra.
%  THCDF     - Marginal wave height, Hd, cdf for Torsethaugen spectra.
%  THVPDF    - Joint (Vcf,Hd) pdf for Torsethaugen spectra.
%
%  NEEDS MORE WORK
%  ~~~~~~~~~~~~~~~~
%  
%  ochi98XXX     - (+) Need a explanation to the scaling.
%  tay81cdf      - (+) Needs testing.
