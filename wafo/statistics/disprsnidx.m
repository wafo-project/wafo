function  [di,b_u,ok_u] = disprsnidx(data,varargin)
%DISPRSNIDX Dispersion Index vs threshold
%
%  CALL [DI,b_u,ok_u] = disprsnidx(data,options)
%
%  DI   = Dispersion index
%  b_u  = threshold where the number of exceedances in a fixed period (Tb) is
%         consistent with a Poisson process.
% ok_u  = all thresholds where the number of exceedances in a fixed period (Tb) is
%         consistent with a Poisson process.
%  data = twocolumn matrix with sampled times in the first column
%               and values the second columns
% options = options structure defining estimation of DI. 
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 10).
%        .umin : Minimum threshold (default umin = sd(1) where sd = sort(data))
%        .umax : Maximum threshold (default umax = sd(end-Nmin) where sd = sort(data); )
%        .Nu   : number of threshold values (default min(length(data),100))
%        .u    : (default linspace(umin,umax,Nu)
%        .Tb   : Block period (same unit as the sampled times)  (default 1)
%        .alpha: Confidence coefficient (default 0.05)
%        .plotflag:
%
%   DISPRSNIDX estimate the Dispersion Index (DI) as function of threshold. 
%   DI measures the homogenity of data and the purpose of DI is to determine 
%   the threshold where the number of exceedances in a fixed period (Tb) is
%   consistent with a Poisson process. For a Poisson process the DI is one. 
%   Thus the threshold should be so high that DI is not significantly
%   different from 1.
%
%   The Poisson hypothesis is not rejected if the estimated DI is between:
%
%   chi2(alpha/2, M-1)/(M-1)< DI < chi^2(1 - alpha/2, M-1 }/(M - 1)
%
%  where M is the total number of fixed periods/blocks -generally
%   the total number of years in the sample.
%
% Example
%  xn = load('sea.dat');
%  Ie = findpot(xn,0,5);
%  [di,u] = disprsnidx(xn(Ie,:),'Tb', 100);
%  plot(di) % a threshold around 1 seems appropriate.
%  vline(u)
%
% See also reslife, fitgenparrange, extremalidx


% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful, but without any warranty; without even
% the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public
% License for more details.
% The GNU General Public License can be obtained from http://www.gnu.org/copyleft/gpl.html. You
% can also obtain it by writing to the Free Software Foundation, Inc., 59 Temple Place – Suite 330, Boston,
% MA 02111-1307, USA.



%References
% Ribatet, M. A.,(2006),
% A User's Guide to the POT Package (Version 1.0)
% month = {August},
% url = {http://cran.r-project.org/}
%
% Cunnane, C. (1979) Note on the poisson assumption in
%     partial duration series model. Water Resource Research, 15\bold{(2)}
%     :489--494.}


% History
% By pab Nov 2007
% translated from diplot in POT package in R by Mathieu Ribatet



options = struct('u',[],'umin',[],'umax',[],'Nu',200,'Nmin',10,'Tb',1,'alpha',0.05,'plotflag',0);

if ischar(data) && strcmpi(data,'defaults')
  di = options;
  return
end
options = parseoptions(options,varargin{:});

ti = data(:,1)-min(data(:,1));

t = floor(ti./options.Tb)+1;
d = (data(:,2));
n = length(d);



if isempty(options.u)
  sd = sort(d);
  if ~isempty(options.Nmin)
    options.Nmin = max(options.Nmin,0);
  end
  if options.Nmin>n/2
    warning('WAFO:DISPRSNIDX','Nmin possibly too large!')
  end
  if (n<=options.Nmin)
    error('Not enough data for a Dispersion Index')
  end

  if isempty(options.umin)
    options.umin = sd(1);
  else
    options.umin = max(options.umin,sd(1));
  end

  sdmax = sd(n-options.Nmin);
  if isempty(options.umax)
    options.umax = sdmax;
  else
    options.umax = min(options.umax,sdmax);
  end

  if isempty(options.Nu)
    options.Nu = min(n-options.Nmin,100);
  end
  u = linspace(options.umin,options.umax,options.Nu).';
else
  u = options.u(:);
end

Nu = length(u);
nan1 = nan;
di1 = nan1(ones(Nu,1));



mint = min(t); % mint should be 0.
maxt = max(t);
M = maxt-mint+1;
occ = zeros(1,M);

for ix = 1:Nu
  excess = d > u(ix);
  lambda = sum(excess) / M;
  for block = 1:M
    occ(block) = sum(excess(t == block));
  end  
  di1(ix) = var(occ)/lambda;
end
p = 1-options.alpha;
diUp = di1;
diLo = di1;
diUp(:) = invchi2(1-options.alpha/2,M-1)/(M-1);
diLo(:) = invchi2(options.alpha/2,M-1)/(M-1);
options.CI = [diLo,diUp];

if nargout>1
  % Find appropriate threshold
  k1 = find(diLo< di1 & di1<diUp);
  if any(k1)
    k = find(mean(di1(k1))<di1(k1),1,'first');
    ok_u = u(k1);
    b_u = ok_u(k);
  else
    b_u = [];
    ok_u = [];
  end
end
  
CItxt = sprintf('%d%s CI',100*p,'%');
titleTxt = 'Dispersion Index plot';

di = createwdata('data',di1,'args',u,...
'dataCI',options.CI,'title',titleTxt,'labels',{'Threshold','Dispersion Index'},...
    'workspace',options,'note',titleTxt,'caption',CItxt);


if ~isoctave
  %di =
  %createpdf('f',[di1,options.CI],'x',{u},'title',titleTxt,'labx',{'Threshold','Dispersion Index'});
  di = wdata(di);
  if options.plotflag>0
    plot(di,options.plotflag,'.')
  end
end

