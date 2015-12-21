function [F,G] = plotedfcnd(z,varargin) 
%PLOTEDFCND Plot Empirical Distribution Function CoNDitioned that X>=c
%          and optionally compare it with distribution G.
%
%  CALL:  [F,Gc] = plotedfcnd(X,c,G,plotflag,sym);
%
%        F  = conditional EDF of X given X>=c, (pdf-struct or wdata).
%        Gc = conditional CDF, G(x;X>c) =  G(x>c)/(1-G(x<c))  (pdf-struct or wdata).
%        X  = data vector.
%        c  = value to be conditioned on (default c = min(x,0)).
%        G  = cdf, two column matrix or FDATA object as 
%            returned from the FITXXX functions(optional).
%           
%  plotflag = 0  no plotting       
%             1 plot cdf F(x|X>=c)   (default)
%               See pdfplot for more details on plotoptions 
%       sym = {s1,s2} cell array or comma separated list of plot 
%             symbols for F and G, respectively.
%             (default {'b','r--'})
% 
% NOTE:  SYM can be given anywhere after X
% 
% Example:
%   x=linspace(0,6,200)';
%   R = rndray(2,100,1);
%   plotedfcnd(R,1,[x,cdfray(x,2)],'g','b') %  plot
%   F = edf(R)
%
% See also edfcnd, plotedf, pdfplot, cumtrapz


% Copyright (C) 2000 WAFO-group
%
% This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Tested on: Matlab 5.3, 5.2, 5.1
% History:
% revised PJ 01-Apr-2001: updated help text
% revised pab 23.06.2000
% - added ih = ishold
% - added sym
% - added varargin
% revised ms 13.06.2000
% - pdf usage removed (can't distinguish between a cdf and an increasing pdf)
% - changed axis labelling in figures
% - changed to staircaseplot when plotflag=2,3
% - revised header info
% revised pab 08.06.2000
% - fixed normalization of f if c>min(z)
% revised pab 07.03.2000
% - enabled so that f may be empty while plotflag is given 
% modified by Per A. Brodtkorb 10.11.98
% to accept both pdf and cdf's. Also enabled new plotting features,
% plotting of  probability of exceedances on a semilogy paper .... 

error(nargchk(1,6,nargin))
ih = ishold;

% default values
%~~~~~~~~~~~~~~~
c        = floor(min(min(z),0));
plotflag = 1; 
G=[];
sym ={[],'r--'}; % default plot symbols for the empirical
                 %  theoretical pdf,respectively
                 
%	  The method (default 3) is
%	  1. Interpolation so that F(X_(k)) == (k-0.5)/n.
%	  2. Interpolation so that F(X_(k)) == k/(n+1).
%	  3. The empirical distribution. F(X_(k)) = k/n
                 
%method = 2;                              
                              

P  = varargin;
Np = length(P);
if Np>0
  strix = zeros(1,Np);
  cellix = strix;
  for ix=1:Np, % finding symbol strings or cell array of symbol strings
    strix(ix)  = ischar(P{ix});
    cellix(ix) = iscell(P{ix});
  end
  k  = find(strix);
  k1 = find(cellix);
  if any(k)
    Nk = length(k);
    if Nk>2,  warning('WAFO:PLOTEDFCND','More than 2 strings are not allowed'),    end
    iy = 1;
    for ix=k      
      sym{iy} = P{ix};
      iy=iy+1;
    end
    Np = Np-Nk;
    P  = {P{find(~strix)}}; % remove strings from input
  elseif any(k1) % cell array of strings
    tmp = P{k1};
    Nk = length(tmp);
    if Nk>2,  warning('WAFO:PLOTEDFCND','More than 2 strings are not allowed'),    end
    for ix=1:min(Nk,2)
      if ~isempty(tmp{ix}) && ischar(tmp{ix}), sym{ix}=tmp{ix};end
    end
    Np = Np-1;
    P  = {P{find(~cellix)}}; % remove cell array of strings from input
  end
  if Np>0 && ~isempty(P{1}), c        = P{1};end
  if Np>1 && ~isempty(P{2}), G        = P{2};end
  if Np>2 && ~isempty(P{3}), plotflag = P{3};end
end

if isempty(sym{1}), 
  sym{1}='b-';
%  if plotflag == 1, sym{1} = 'b-'; else sym{1}='b.'; end
end

[Fz,G] = edfcnd(z,c,G,'wdata',true);

if plotflag>0
  plotflag1 = plotflag - mod(plotflag,10)+2;
  plot(Fz,plotflag1,sym{1});
  if ~isempty(G)
    hold on, plot(G,plotflag,sym{2});
  end
  if ~ih, hold off,end
end
if nargout>0
  F = Fz;
end

return
