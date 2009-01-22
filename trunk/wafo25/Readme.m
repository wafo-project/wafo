%README New features, bug fixes, and changes in WAFO Toolbox - Wave Analysis for Fatigue and Oceanography 
% Version 2.5  08-Jan-2009
%
%  NEW DIRECTORIES
% ~~~~~~~~~~~~~~~~
%  GRAPHUTIL : directory with many useful plott function
% NUMDIFFTOOLS : directory with tools for automatic numerical differentiation of functions
% POLYUTIL : directory tools for evaluating polynomials.
%
%  NEW FUNCTIONS 	
%  ~~~~~~~~~~~~~~
% Spectrum estimation:
%   welch_psd, welch_psd2, welch_cpsd: translated from Octave (GPL)
%
%  New distribution functions added:
%  Poisson, cauchy, hypergeometric, 
%
% Profile log likelihood functions
%   lnkexp            - Link for x,F and parameters of Exponential distribution
%   lnkgenpar         - Link for x,F and parameters of Generalized Pareto distribution
%   lnkgev            - Link for x,F and parameters of Generalized Extreme value distribution
%   lnkgumb           - Link for x,F and parameters of Gumbel distribution
%   lnkgumbtrnc       - Link for x,F and parameters of truncated Gumbel distribution
%   lnkray            - Link for x,F and parameters of Rayleigh distribution
%   lnkweib           - Link for x,F and parameters of the Weibull distribution
%   logps             - Moran's negative log Product Spacings statistic
%   ciproflog         - Confidence Interval using Profile Log- likelihood or Product Spacing- function.
%   proflog           - Profile Log- likelihood or Product Spacing-function.
%   FINDCIPROFLOG     - Find Confidence Interval from proflog function
%
%   decluster         - Decluster peaks over threshold values
%   extremalidx       - Extremal Index measuring the dependence of data
%   findpot           - Find indices to Peaks over threshold values
%   prb2retper        - Return period from Probability of exceedance. 
%   retper2prb        - Probability of exceedance from return period.
%
%   fitgenparrange    - Parameter estimates for GPD model vs thresholds
%   disprsnidx        - Dispersion Index vs threshold
%   reslife           - Mean Residual Life, i.e., mean excesses vs thresholds
%   plotdisprsnidx    - Plot Dispersion Index vs thresholds
%   plotreslife       - Plot Mean Residual Life (mean excess vs thresholds)
%  
%   logit             - Logit function.
%   logitinv          - Inverse logit function.
%   regglm            - Generalized Linear Model regression
%   reglm             - Fit multiple Linear Regression Model.
%   reglogit          - Fit ordinal logistic regression model.
%   regnonlm          - Non-Linear Model Regression  
%   regsteplm         - Stepwise predictor subset selection for Linear Model regression
%
%   princomp          -  Compute principal components of X
%
%   ranktrf           - Rank transformation of data material.
%   spearman          - Spearman's rank correlation coefficient.
%   var2corr          - Variance matrix to correlation matrix conversion.
%   lmoment           - L-moment based on order statistics.
%
%   plotbox           - Plot box-and-whisker diagram
%   plotdensity       - Plot density.
%   plotfitsumry      - Plot diagnostic of fit to data
%
%   anovan            - multi-way analysis of variance (ANOVA)
%   testgumb          - Tests if shape parameter in a GEV is equal to zero
%   testmean1boot     - Bootstrap t-test for the mean equal to 0.
%   testmean1n        - Test for mean equals 0 using one-sample T-test
%   testmean2n        - Two-sample t-test for mean(x) equals mean(y)
%   testmean1r        - Wilcoxon signed rank test for H0: mean(x) equals 0.
%   testmean2r        - Wilcoxon rank-sum test for H0: mean(x) equals mean(y).
%
% Confidence interval estimation
%   ciboot            - Bootstrap confidence interval.
%   ciquant           - Nonparametric confidence interval for quantile
%   momci1b           - Moment confidence intervals using Bootstrap
%
% Bootstrap & jacknife estimates
%   covboot           - Bootstrap estimate of the variance of a parameter estimate.
%   covjack           - Jackknife estimate of the variance of a parameter estimate.
%   stdboot           - Bootstrap estimate of the parameter standard deviation.
%   stdjack           - Jackknife estimate of the standard deviation of a parameter 
%  
%   pdfdiscrete       - Discrete PDF
%   pdfempirical      - Empirical pdf
%   cdfdiscrete       - Discrete cdf
%   cdfempirical      - Empirical cdf
%   invdiscrete       - Disrete quantile
%   invempirical      - Empirical quantile. 
%  
%
%  CHANGES
%  ~~~~~~~~~
%  XXXpdf, XXXcdf, XXXrnd, XXXinv, XXXfit, XXXplot renamed to pdfXXX, cdfXXX, rndXXX,
%  invXXX fitXXX, plotXXX, respectively.
%  XXXstat renamed to momXXX
%  fitXXX now returns a structure instead of a vector.
%  All pdfXXX, cdfXXX, invXXX, momXXX now accepts a structure as returned
%  from fitXXX functions.
%  All momXXX functions return mean, std, skewness and kurtosis.
%  iqr      renamed to iqrange
%  identify renamed to clickslct
%  whisto   renamed to histgrm
%  mle      renamed to mlest
%  empdistr renamed to plotedf
%  cempdistr renamed to plotedfcnd
%  wskewness renamed to skew 
%  wkurtosis renamed to kurt 
%  wquantile renamed to percentile
%
% Version 2.1.1   07-Sep-2005
%
%  FIXES
%  ~~~~~
%   o Some minor bugfixes in the fatigue part. 
%   o Removed all ':' from the 'See also' lines in all m-files, because
%     the ':' obstructed the possibility for clicking on the links in the
%     'See also' line in matlab 7. 
%   o Many new wavemodels are available in the WAVEMODELS directory.
%   o Some minor bugfixes in WSTATS, SPEC, KDETOOLS DATA and MISC
%     directories.
%   o The HTML based documentation can be launched from matlab by typing
%     wafohelp.
%
%   NEW FUNCTIONS
%   ~~~~~~~~~~~~~
%   o New functions for computing level crossings for noncentral Chi^2
%     processes in TRGAUSS directory.
%   o Some new functions in the CYCLES and WSIM directories
%
%-------------------------------------------------------------------
% Version 2.1.0   07-Apr-2004
%
% WAFO has been upgraded. Major changes concern Chapter 3 in the 
% tutorial. The calls there are not valid anymore. 
% 
% A new HTML based documentation is available, with for example
% cross references between calls of functions. 
% 
% For PC Windows, the executable programs have been replaced by 
% mex files (these are not compiled for other systems). 
%--------------------------------------------------------------------
% Version 2.0.5   14-Jun-2003
%
% This version contains all the reported bugs and updates on the WAFO web.
% It also contains additional bug fixes, and new/updated routines.
% 
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
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
%  FIXES
%  ~~~~~
%  See the Readme files in the modules.
%
%  ENHANCEMENTS
%  ~~~~~~~~~~~~
%  See the Readme files in the modules.
% 
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%  See the Readme files in the modules.
%  
%  SHORT HISTORY
%  ~~~~~~~~~~~~~
%  The base for the WAFO toolbox was the Wave Analysis Toolbox (WAT)
%  and the Fatigue Analysis Toolbox (FAT). 
%  Several new functions exist along with enhancements of the old 
%  WAT and FAT functions. The enhancements also include a more consistent 
%  naming convention for the m-files.
%  Most of the functions in WAFO are incompatible with the functions
%  of WAT and FAT.
