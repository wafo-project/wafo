function y = stirlerr(n)
%STIRLERR Computes  log(n!) - log( sqrt(2*pi*n)*(n/exp(1))^n )
%
%  CALL stirlerr(n)
%
% STIRLERR computes the error of the Stirling approximation, i.e.,
%  log(n!) - log( sqrt(2*pi*n)*(n/exp(1))^n
%
% Example
%  assert(stirlerr(1), 0.0810614667953273, 1e-12)
%  assert(stirlerr(20), 0.00416631969199693, 1e-13)
%  assert(stirlerr(60), 0.00138887602982701, 1e-15)
%  assert(stirlerr(100), 8.33330555634921e-004, 1e-16)
%  assert(stirlerr(1000), 8.33333305555556e-005, 1e-17)
%  assert(stirlerr(realmax), 0, 1e-17)
%
% See also binom, pdfbin, pdfpois


% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.
% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }


S0 = 0.083333333333333333333; % /* 1/12 */
S1 = 0.00277777777777777777778;% /* 1/360 */
S2 = 0.00079365079365079365079365;% /* 1/1260 */
S3 = 0.000595238095238095238095238;% /* 1/1680 */
S4 = 0.0008417508417508417508417508;% /* 1/1188 */

% sfe=[0, 0.081061466795327258219670264,...
% 0.041340695955409294093822081, 0.0276779256849983391487892927,...
% 0.020790672103765093111522771, 0.0166446911898211921631948653,...
% 0.013876128823070747998745727, 0.0118967099458917700950557241,...
% 0.010411265261972096497478567, 0.0092554621827127329177286366,...
% 0.008330563433362871256469318, 0.0075736754879518407949720242,...
% 0.006942840107209529865664152, 0.0064089941880042070684396310,...
% 0.005951370112758847735624416, 0.0055547335519628013710386899];

%y = n;
%nL16    = n<16;
%y(nL16) = sfe(floor(n(nL16))+1);
n = min(n, realmax*1e-3);
y = gammaln(n+1) - log( sqrt(2*pi.*n).*(n./exp(1)).^n );


nn = double(n);
nn = nn.*nn;

n500    = 500<n;
y(n500) = (S0-S1./nn(n500))./n(n500);
n80     = 80<n & n<=500;
y(n80)  = (S0-(S1-S2./nn(n80))./nn(n80))./n(n80);
n35     = 35< n & n<=80;
if any(n35)
  nn35   = nn(n35);
  y(n35) = (S0-(S1-(S2-S3./nn35)./nn35)./nn35)./n(n35);
end
n15      = 15< n & n<=35;
if any(n15)
  nn15   = nn(n15);
  y(n15) = (S0-(S1-(S2-(S3-S4./nn15)./nn15)./nn15)./nn15)./n(n15);
end

%!test assert(stirlerr(1), 0.0810614667953273, 1e-12)
%!test assert(stirlerr(20), 0.00416631969199693, 1e-13)
%!test assert(stirlerr(60), 0.00138887602982701, 1e-15)
%!test assert(stirlerr(100), 8.33330555634921e-004, 1e-16)
%!test assert(stirlerr(1000), 8.33333305555556e-005, 1e-17)
%!test assert(stirlerr(realmax), 0, 1e-17)