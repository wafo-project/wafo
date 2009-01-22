function [h,ep] = mk87rnd(Hrms,Erms, N)
%MK87RND Random points from MK87 distribution of steepness and wave height.
%
% CALL:  [H E] = mk87rnd(Hrms,Erms,N )
%
%  H,E        = two vectors of simulated data of length N
%  Hrms, Erms =  parameters of the MK87 distribution
%
% Example:
%  Hs = 7;Tz=10;  
%  Hrms = 0.715*Hs;
%  Erms = 0.0202+0.826*Hs/(Tz^2);
%  h = linspace(0,3*Hrms)'; s = linspace(0,5*Erms)';
%  f = mk87pdf2(h,s,Hrms,Erms);
%  pdfplot(f)
%  [H E] = mk87rnd(Hrms,Erms,500);
%  hold on, plot(E,H,'.'), hold off
%
% See also  mk87pdf, rndweib, rndlog

%   References:
%   Myrhaug, D. and Kjelsen S.P. (1987) 
%  'Prediction of occurences of steep and high waves in deep water'.
%   Journal of waterway, Port, Coastal and Ocean Engineers, Vol. 113, pp 122--138
%
%   Myrhaug & Dahle (1984) Parametric modelling of joint probability
%   density distributions for steepness and asymmetry in deep water 
%   waves

% tested on: matlab 5.1
% history:
% revised pab 01.04.2001
% -added example
% revised pab 04.11.2000
% no dependence on stats toolbox anymore
% by  Per A. Brodtkorb 19.11.1998



if nargin <  2 
  error('Requires at least two input arguments.'); 
end
  

if nargin < 2||isempty(Erms),
  Erms=1;
end
if nargin < 1||isempty(Hrms),
  Hrms=1;
end
[icode Erms,Hrms] = iscomnsize(Erms,Hrms);
if ~icode
  error('Requires non-scalar arguments to match in size.');
end
[n m]=size(Hrms);
Hrms=Hrms(:);Erms=Erms(:);
if nargin < 3,
  N=n*m;
end
%h=zeros(N,1);
%ep=h;
  


% NB! weibpdf must be modified to correspond to
% pdf=x^(b-1)/a^b*exp(-(x/a)^b) or else insert
% weibpdf=2.39.*h.^1.39/(1.05^2.39).*exp(-(h./1.05).^2.39);

h=rndweib(1.05,2.39,N,1);
sig=(-0.21*atan(2*(h-1.4))+0.325);
ep=rndlognorm(my(h),sig).*Erms;
h=h.*Hrms;

function y=my(h)
y=zeros(size(h));
ind=(h <= 1.7);
h1=h(ind);
y(ind)=0.024-1.065.*h1+0.585.*h1.^2;
%ind=find(h > 1.7);
%h2=h(~ind);
y(~ind)=0.32*atan(3.14*(h(~ind)-1.7))-0.096;
return
