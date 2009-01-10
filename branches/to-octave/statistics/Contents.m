% Module STATISTICS in WAFO Toolbox. 
% Version 2.2.1  07-Dec-2007
%
% What's new
%   Readme            - New features, bug fixes, and changes in STATISTICS. 
%
% Parameter estimation
%   fitbeta           - Parameter estimates for Beta data.
%   fitchi2           - Parameter estimates for Chi squared data.
%   fitexp            - Parameter estimates for Exponential data.
%   fitgam            - Parameter estimates for Gamma data.
%   fitgengam         - Parameter estimates for Generalized Gamma data.
%   fitgenpar         - Parameter estimates for Generalized Pareto data
%   fitgenparml       - Internal routine for fitgenpar (ML estimates for GPD data) 
%   fitgenparrange    - Parameter estimates for GPD model over a range of thresholds
%   fitgev            - Parameter estimates for GEV data
%   fitgumb           - Parameter estimates for Gumbel data.
%   fitinvnorm        - Parameter estimates for Inverse Gaussian data.
%   fitlognorm        - Parameter estimates for Lognormal data.
%   fitmarg2d         - Parameter estimates for MARG2D data.
%   fitmargcnd2d      - Parameter estimates for DIST2D data.
%   fitnorm           - Parameter estimates for Normal data.
%   fitray            - Parameter estimates for Rayleigh data.
%   fitraymod         - Parameter estimates for Truncated Rayleigh data.
%   fitt              - Parameter estimates for Student's T data.
%   fitweib           - Parameter estimates for Weibull data.
%   fitweib2d         - Parameter estimates for 2D Weibull data.
%   fitweibmod        - Parameter estimates for truncated Weibull data.
%
%   likgenpar         - Log likelihood function for GPD
%   likweib2d         - 2D Weibull log-likelihood function.
%
%   loglike           - Negative Log-likelihood function.
%   logps             - Moran's negative log Product Spacings statistic
%   mlest             - Maximum Likelihood or Maximum Product Spacing estimator
%
% Probability density functions (pdf)
%   pdfbeta           - Beta probability density function
%   pdfbin            - Binomial probability mass function
%   pdfcauchy         - Cauchy's probability density function
%   pdfchi2           - Chi squared probability density function
%   pdfdiscrete       - Discrete PDF
%   pdfempirical      - Empirical pdf
%   pdfexp            - Exponential probability density function
%   pdff              - Snedecor's F probability density function
%   pdffrech          - Frechet probability density function
%   pdfgam            - Gamma probability density function
%   pdfgengam         - Generalized Gamma probability density function
%   pdfgengammod      - Modified Generalized Gamma probability density function (stable)
%   pdfgenpar         - Generalized Pareto probability density function
%   pdfgev            - Generalized Extreme Value probability density function
%   pdfgumb           - Gumbel probability density function.
%   pdfhyge           - Hypergeometric probability mass function
%   pdfinvnorm        - Inverse Gaussian probability density function
%   pdflognorm        - Lognormal probability density function
%   pdfmarg2d         - Joint 2D PDF due to Plackett given as  f{x1}*f{x2}*G(x1,x2;Psi). 
%   pdfmargcnd2d      - Joint 2D PDF computed as f(x1|X2=x2)*f(x2)
%   pdfnorm           - Normal probability density function
%   pdfnorm2d         - Bivariate Gaussian distribution  
%   pdfnormnd         - Multivariate Normal probability density function.  
%   pdfray            - Rayleigh probability density function
%   pdfraymod         - Truncated Rayleigh probability density function
%   pdft              - Student's T probability density function
%   pdfpois           - Poisson probability mass function
%   pdfweib           - Weibull probability density function
%   pdfweib2d         - 2D Weibull probability density function (pdf).
%   pdfweibmod        - Truncated Weibull probability density function
%
% Cumulative distribution functions (cdf)
%   cdfcauchy         - Cauchy cumulative distribution function
%   cdfdiscrete       - Discrete cdf
%   cdfempirical      - Empirical cdf
%   cdfmarg2d         - Joint 2D CDF due to Plackett 
%   cdfmargcnd2d      - Joint 2D CDF computed as int F(X1<v|X2=x2).*f(x2)dx2
%   cdfmargcnd2dfun   - is an internal function to cdfmargcnd2d prbmargcnd2d. 
%   cdfnormnd         - Multivariate normal cumulative distribution function (cdf).
%   cdfweib2d         - Joint 2D Weibull cumulative distribution function  
%   cdfbeta           - Beta cumulative distribution function
%   cdfbin            - Binomial cumulative probability function
%   cdfchi2           - Chi squared cumulative distribution function
%   cdfexp            - Exponential cumulative distribution function
%   cdff              - Snedecor's F cumulative distribution function
%   cdffrech          - Frechet cumulative distribution function
%   cdfgam            - Gamma cumulative distribution function
%   cdfgengam         - Generalized Gamma cumulative distribution function
%   cdfgengammod      - Modified Generalized Gamma cumulative distribution function
%   cdfgenpar         - Generalized Pareto cumulative distribution function
%   cdfgev            - Generalized Extreme Value cumulative distribution function
%   cdfgumb           - Gumbel cumulative distribution function.
%   cdfhyge           - The hypergeometric cumulative probability function
%   cdfinvnorm        - Inverse Gaussian cumulative distribution function
%   cdflognorm        - Lognormal cumulative distribution function
%   cdfmargcnd2d      - Joint 2D CDF computed as int F(X1<v|X2=x2).*f(x2)dx2
%   cdfnorm           - Normal cumulative distribution function
%   cdfray            - Rayleigh cumulative distribution function
%   cdfraymod         - Modified Rayleigh cumulative distribution function
%   cdft              - Student's T  cumulative distribution function
%   cdfpois           - Poisson cumulative distribution function
%   cdfweib           - Weibull cumulative distribution function
%   cdfweibmod        - Truncated Weibull cumulative distribution function
%
%   edf               - Empirical Distribution Function 
%   edfcnd            - Empirical Distribution Function CoNDitioned that X>=c
%
%   prbmargcnd2d      - returns the probability for rectangular regions.
%   prbweib2d         - returns the probability for rectangular regions.
%   margcnd2dsmfun    - Smooths the MARGCND2D distribution parameters. 
%   margcnd2dsmfun2   - Smooths the MARGCND2D distribution parameters. 
%
% Inverse cumulative distribution functions
%   invbeta           - Inverse of the Beta distribution function
%   invbin            - Inverse of the Binomial distribution function.
%   invcauchy         - Inverse of the Cauchy distribution function.
%   invchi2           - Inverse of the Chi squared distribution function
%   invcmarg2d        - Inverse of the conditional cdf of X2 given X1.
%   invcweib2d        - Inverse of the conditional 2D weibull cdf of X2 given X1.
%   invdiscrete       - Disrete quantile
%   invempirical      - Empirical quantile. 
%   invexp            - Inverse of the Exponential distribution function
%   invf              - Inverse of the Snedecor's F distribution function
%   invfrech          - Inverse of the Frechet distribution function
%   invgam            - Inverse of the Gamma distribution function
%   invgengam         - Inverse of the Generalized Gamma distribution function
%   invgengammod      - Inverse of the Generalized Gamma distribution function
%   invgenpar         - Inverse of the Generalized Pareto distribution function
%   invgev            - Inverse of the Generalized Extreme Value distribution function
%   invgumb           - Inverse of the Gumbel distribution function.
%   invhyge           - Inverse of the Hypergeometric distribution function.
%   invinvnorm        - Inverse of the Inverse Gaussian distribution function
%   invlognorm        - Inverse of the Lognormal distribution function
%   invnorm           - Inverse of the Normal distribution function
%   invray            - Inverse of the Rayleigh distribution function
%   invt              - Inverse of the Student's T distribution function
%   invweib           - Inverse of the Weibull distribution function
%   invpois           - Inverse of the Poisson distribution function
%   invraymod         - Inverse of the modified Rayleigh distribution function
%   invweibmod        - Inverse of the modified Weibull distribution function
%
% Random number generators
%   rndalpha          - Random matrices from a symmetric alpha-stable distribution
%   rndbeta           - Random matrices from a Beta distribution
%   rndbin            - Random numbers from the binomial distribution
%   rndboot           - Simulate a bootstrap resample from a sample.
%   rndcauchy         - Random matrices a the Cauchy distribution.
%   rndchi2           - Random matrices from a Chi squared distribution.
%   rnddiscrete       - Random sample
%   rndempirical      - Bootstrap sample
%   rndexp            - Random matrices from an Exponential distribution
%   rndf              - Random matrices from the Snedecor's F distribution
%   rndfrech          - Random matrices from a Frechet distribution.
%   rndgam            - Random matrices from a Gamma distribution.
%   rndgengam         - Random matrices from a Generalized Gamma distribution.
%   rndgengammod      - Random matrices from a Generalized Modified Gamma distribution.
%   rndgenpar         - Random matrices from a Generalized Pareto Distribution
%   rndgev            - Random matrices from a Generalized Extreme Value distribution
%   rndgumb           - Random matrices from a Gumbel distribution.
%   rndhyge           - Random numbers from the Hypergeometric distribution
%   rndinvnorm        - Random matrices from a Inverse Gaussian distribution.
%   rndlognorm        - Random matrices from a Lognormal distribution.
%   rndmarg2d         - Random points from a MARG2D distribution 
%   rndmargcnd2d      - Random points from a MARGCND2D distribution 
%   rndnorm           - Random matrices from a Normal distribution.
%   rndnormnd         - Random vectors from a multivariate Normal distribution
%   rndpois           - Random matrices from a Poisson distribution
%   rndray            - Random matrices from a Rayleigh distribution
%   rndraymod         - Random matrices from modified Rayleigh distribution
%   rndt              - Random matrices from a Student's T distribution
%   rndweib           - Random matrices a the Weibull distribution.
%   rndweibmod        - Random matrices from the modified Weibull distribution.
%   rndweib2d         - Random numbers from the 2D Weibull distribution.
%
% Moments
%   mombeta           - Mean and variance for the Beta distribution.
%   mombin            - Mean and variance for the BINOMIAL distribution.
%   momchi2           - Mean and variance for the Chi squared distribution.
%   momexp            - Mean and variance for the Exponential distribution.
%   momf              - Mean and variance for the Snedecor's F distribution.
%   momfrech          - Mean and variance for the Frechet distribution.
%   momgam            - Mean and variance for the Gamma distribution.
%   momgengam         - Mean and variance for the Generalized Gamma distribution.
%   momgenpar         - Mean and variance for the Generalized Pareto distribution.
%   momgev            - Mean and variance for the GEV distribution.
%   momgumb           - Mean and variance for the Gumbel distribution.
%   momhyge           - Mean and variance for the Hypergeometric distribution.
%   mominvnorm        - Mean and variance for the Inverse Gaussian distribution.
%   momlognorm        - Mean and variance for the Lognormal distribution.
%   mommarg2d         - Mean and variance for the MARG2D distribution.
%   mommargcnd2d      - Mean and variance for the MARGCND2D distribution
%   momnorm           - Mean and variance for the Normal distribution.
%   mompois           - Mean and variance for the Poisson distribution.
%   momray            - Mean and variance for the Rayleigh distribution.
%   momt              - Mean and variance for the Student's T  distribution.
%   momweib           - Mean and variance for the Weibull  distribution.
%   momweib2d         - Mean and variance for the 2D Weibull distribution 
%
% Profile log likelihood functions
%   lnkexp            - Link for x,F and parameters of Exponential distribution
%   lnkgenpar         - Link for x,F and parameters of Generalized Pareto distribution
%   lnkgev            - Link for x,F and parameters of Generalized Extreme value distribution
%   lnkgumb           - Link for x,F and parameters of Gumbel distribution
%   lnkgumbtrnc       - Link for x,F and parameters of truncated Gumbel distribution
%   lnkray            - Link for x,F and parameters of Rayleigh distribution
%   lnkweib           - Link for x,F and parameters of the Weibull distribution
%   loglike           - Negative Log-likelihood function.
%   logps             - Moran's negative log Product Spacings statistic
%   ciproflog         - Confidence Interval using Profile Log- likelihood or Product Spacing- function.
%   proflog           - Profile Log- likelihood or Product Spacing-function.
%   findciproflog     - Find Confidence Interval from proflog function
%
% Extremes 
%   decluster         - Decluster peaks over threshold values
%   extremalidx       - Extremal Index measuring the dependence of data
%   findpot           - Find indices to Peaks over threshold values
%   fitgev            - Parameter estimates for GEV data
%   fitgenpar         - Parameter estimates for Generalized Pareto data
%   prb2retper        - Return period from Probability of exceedance.  
%   retper2prb        - Probability of exceedance from return period. 
%
% Threshold selection
%   fitgenparrange    - Parameter estimates for GPD model vs thresholds
%   disprsnidx        - Dispersion Index vs threshold
%   reslife           - Mean Residual Life, i.e., mean excesses vs thresholds
%   plotdisprsnidx    - Plot Dispersion Index vs thresholds
%   plotreslife       - Plot Mean Residual Life (mean excess vs thresholds)
%   
% Regression models
%   logit             - Logit function.
%   logitinv          - Inverse logit function.
%   regglm            - Generalized Linear Model regression
%   reglm             - Fit multiple Linear Regression Model.
%   reglogit          - Fit ordinal logistic regression model.
%   regnonlm          - Non-Linear Model Regression   
%   regsteplm         - Stepwise predictor subset selection for Linear Model regression
%
% Factor analysis
%   princomp          -  Compute principal components of X
%
% Descriptive Statistics
%   ranktrf           - Rank transformation of data material.
%   spearman          - Spearman's rank correlation coefficient.
%   mean              - Computes sample mean (in matlab toolbox)
%   median            - Computes sample median value (in matlab toolbox)
%   std               - Computes standard deviation (in matlab toolbox)
%   var               - Computes sample variance (in matlab toolbox)
%   var2corr          - Variance matrix to correlation matrix conversion.
%   cov               - Computes sample covariance matrix (in matlab toolbox)
%   corrcoef          - Computes sample correlation coefficients (in matlab toolbox)
%   skew              - Computes sample skewness
%   kurt              - Computes sample kurtosis
%   lmoment           - L-moment based on order statistics.
%   percentile        - Empirical quantile (percentile).
%   iqrange           - Computes the Inter Quartile Range
%   range             - Computes the range between the maximum and minimum values. 
%
% Statistical plotting
%   clickslct         - Select points in a plot by clicking with the mouse.
%   histgrm           - Plot histogram
%   plotbox           - Plot box-and-whisker diagram
%   plotdensity       - Plot density.
%   plotexp           - Plot data on Exponential distribution paper
%   plotedf           - Plot Empirical Distribution Function  
%   plotedfcnd        - Plot Empirical Distribution Function CoNDitioned that X>=c
%   plotfitsumry      - Plot diagnostic of fit to data
%   plotgumb          - Plot data on Gumbel distribution paper.
%   plotkde           - Plot kernel density estimate of PDF
%   plotmarg2dcdf     - Plot conditional CDF of X1 given X2=x2
%   plotmarg2dmom     - Plot conditional mean and standard deviation.
%   plotmargcnd2dcdf  - Plot conditional empirical CDF of X1 given X2=x2
%   plotmargcnd2dfit  - Plot parameters of the conditional distribution
%   plotmargcnd2dmom  - Plot conditional mean and standard deviation
%   plotnorm          - Plot data on a Normal distribution paper
%   plotqq            - Plot empirical quantile of X vs empirical quantile of Y
%   plotray           - Plot data on a Rayleigh distribution paper
%   plotresprb        - Plot Residual Probability.
%   plotresq          - Plot Residual Quantile.
%   plotscatr         - Pairwise scatter plots.
%   plotweib          - Plot data on a Weibull distribution paper
%   plotweib2dcdf     - Plot conditional empirical CDF of X1 given X2=x2
%   plotweib2dmom     - Plot conditional mean and standard deviation
%
% Hypothesis Tests
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
%   stdboot           - Bootstrap estimate of the standard deviation of a parameter.
%   stdjack           - Jackknife estimate of the standard deviation of a parameter 
%
% Design of Experiments
%   yates             - Calculates main and interaction effects using Yates' algorithm.
%   ryates            - Reverse Yates' algorithm to give estimated responses
%   fitmodel          - Fits response by polynomial
%   alias             - Alias structure of a fractional design.
%   cdr               - Complete Defining Relation
%   cl2cnr            - Column Label to Column Number
%   cnr2cl            - Column Number to Column Label.
%   ffd               - Two-level Fractional Factorial Design
%   getmodel          - Return the model parameters
%   sudg              - Some Useful Design Generators
%   plotresponse      - Cubic plot of responses
%   nplot             - Normal probability plot of effects
%
% Misc
%   comnsize          - Calculates common size of all non-scalar arguments.
%   dgammainc         - Incomplete gamma function with derivatives.
%   gammaincln        - Logarithm of incomplete gamma function.
%   parsestatsinput   - Parses inputs to pdfxx, prbxx, invxx and rndxx functions.
%   createfdata       - Distribution parameter struct constructor
%   getdistname       - Return the distribution name
%
%   stdize            - Standardize columns to have mean 0 and standard deviation 1.
%   center            - Recenter columns to have mean 0.
%
%  Demo
%   demofitgenpar     - Script to check the variance of estimated parameters






