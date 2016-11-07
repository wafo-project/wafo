function R = dat2cov(xn,varargin)
%DAT2COV Estimate auto covariance function from data.
%
% CALL:  R = dat2cov(x,L,plotflag,dT,flag)
%
%         R = a structure containing:
%                 R = ACF vector length L+1  
%                 t = time lags  length L+1  
%             stdev = estimated large lag standard deviation of the estimate
%                     assuming x is a Gaussian process:
%                     if R(k)=0 for all lags k>q then an approximation 
%                     of the variance for large samples due to Bartlett
%                     var(R(k))=1/N*(R(0)^2+2*R(1)^2+2*R(2)^2+ ..+2*R(q)^2)
%                     for  k>q and where  N=length(x). Special case is
%                     white noise where it equals R(0)^2/N for k>0
%                 h = water depth (default inf)
%                tr = [], transformation 
%              type = 'none'
%              norm = 0 indicating that R is not normalized
%
%        x  = a column data vector or
%             two column data matrix with sampled times and values.
%        L  = the maximum time-lag for which the ACF is estimated.
%	           (Default L=n-1)        
% plotflag = 1 then the ACF is plotted vs lag
%            2 then the ACF is plotted vs lag in seconds
%            3 then the ACF is plotted vs lag and vs lag (sec
%       dT = time-step between data points (default xn(2,1)-xn(1,1) or 1 Hz).
%     flag = 'biased'  : scales the raw cross-correlation by 1/n. (default)
%            'unbiased': scales the raw correlation by 1/(n-abs(k)), 
%                        where k is the index into the result.
%
% Note:  - flag may be given anywhere after x.
%        - x may contain NaN's  (i.e. missing values).
%
% Example: 
%   x = load('sea.dat');
%   rf = dat2cov(x,150,2);
%   assert(rf.R(1:3)', ...
%     [ 0.223686369437041, 0.208384728063068, 0.171107334682617], 1e-10);
%
%   close all;
%
% See also  dat2cor, covplot

%(If you have got the signal toolbox you may modify the dat2cor file)%--------

% Tested on: Matlab 6.0, 5.3, 5.2, 5.1
% History:
% revised jr 02.04.2001
%  - added example, updating help 
% revised pab 09.10.2000
%  - added dat2corchk
%  - updated documentation
%  - added flag option
% revised pab 12.08.99
%  - updated to handle covariance class/structure
% modified by Per A. Brodtkorb 28.10.98 
%  removed some code and replaced it wih a call to dat2cor


[x,L,plotflag,dT,flag,n] = dat2corchk(varargin,xn,nargout);

 
% ind     = find(~isnan(x));
%n0      = length(ind);
%c0      = (n0-1)/n0*(std(x(ind)).^2);% allowing missing data
[R, c0] = dat2cor(x,L,0,dT,flag);
R.R     = R.R*c0;
R.stdev = R.stdev*c0;
R.norm  = 0;

if ((plotflag>0) ),
  covplot(R,L,plotflag);
end


function [x,L,plotflag,dT,flag,n] = dat2corchk(P,xn,Na)
% DAT2CORCHK Helper function for dat2cor.
%
% CALL  [x2,L,plotflag,dT,flag,n]=dat2corchk(P,x1) 
%
%   P = the cell array P of input arguments (between 0 and 7 elements)
%   x = a vector.


x=xn;
[n m]= size(x);

if n<m
 b=m;m=n;n=b; 
 x=x';
end

if n<2, 
  error('The vector must have more than 2 elements!');
end

% Default values
%~~~~~~~~~~~~~~~~  
flag = 'biased';
plotflag = 0;
dT=[]; % sampling frequency unknown
L=n-1;
   % if n<400
%    L=floor(n/2);
%  else
%    L=floor(12*sqrt(n));
%  end

Np=length(P);
strix=zeros(1,Np);
for ix=1:Np, % finding symbol strings 
 strix(ix)=ischar(P{ix});
end


k=find(strix);
if any(k)
  flag = P{k(1)};
  Np   = Np-1;
  P={P{find(~strix)}};
end

if Np>0 && ~isempty(P{1}), L=P{1}; end
if Np>1 && ~isempty(P{2}), plotflag=P{2}; end
if Np>2 && ~isempty(P{3}), dT=P{3}; end

if (Na == 0)&&plotflag==0
  plotflag=3;
end
if plotflag>3, 
  error('Invalid option. Plotflag can only be 0,1,2 or 3');
end
switch m,
 case 1, x=x(:); if isempty(dT),dT=1;end
 case 2, dT=x(2,1)-x(1,1);x=x(:,2);
 otherwise, 
    error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ');
end
return
