function [F,G] = edfcnd(z,c,G,varargin) 
%EDFCND Empirical Distribution Function CoNDitioned that X>=c
%          and optionally computes it for distribution G.
%
%  CALL:  [Fc,Gc] = edfcnd(X,c,G,options);
%
%        Fc = conditional EDF of X given X>=c, (2 column matrix or wdata).
%        Gc = conditional CDF, G(x;X>c) =  G(x>c)/(1-G(x<c))  (2 column matrix or wdata).
%        X  = data vector.
%        c  = value to be conditioned on (default c = min(x,0)).
%        G  = cdf, two column matrix or FDATA object as 
%            returned from the FITXXX functions(optional).
%   options = options struct with fieldnames
%     .method : The method  is
%	            1. Interpolation so that F(X_(k)) == (k-0.5)/n.
%	            2. Interpolation so that F(X_(k)) == k/(n+1).    (default)
%	            3. The empirical distribution. F(X_(k)) = k/n
%     .wdata : If TRUE return as wdata object 
%              otherwise as two column matrix (default)
% 
% Example:
%   x=linspace(0,6,200)';
%   R = rndray(2,100,1);
%   [Fc,Gc] = edfcnd(R,1,[x,cdfray(x,2)]) %  plot
%   plot(Fc),hold on, plot(Gc,'r'), hold off
%   F = edf(R)
%
% See also edf, pdfplot, cumtrapz


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
% by pab2007
% based on old empdistr.

options = struct('method',2,'wdata',false);
error(nargchk(1,inf,nargin))
options = parseoptions(options,varargin{:});

% default values
%~~~~~~~~~~~~~~~

if nargin<2 || isempty(c)
  c = floor(min(min(z),0));
end                
if nargin<3 
  G = [];
end
                      


z = sort(z);
I = find(z>=c);
if isempty(I),  
  error('No data points  z  with  z>=c.'), 
end

z = z(I);
N = length(z);

switch options.method
  case 1
    Fz1 = (0.5:N-0.5)'/N;
  case 3
    Fz1 = (1:N)'/N;
  otherwise
    Fz1 = (1:N)'/(N+1);
end


if c~=-inf,
  ylab = ['F(x| X>=' num2str(c) ')'];
else
  ylab = 'F(x)';
end
xlab = 'x';


if  ~isempty(G)
  if isnumeric(G)
    I = find(G(:,1)>=c);
    if isempty(I),
      error('The cdf  must be defined for at least one value >=c.'),
    end
  
    i = min(I);
    if i > 1 && c>-inf,% Normalize the CDF
      fc = G(i-1,2)+(G(i,2)-G(i-1,2))/(G(i,1)-G(i-1,1))*(c-G(i-1,1));
      G = [c G(I,1)' ; 0  ( G(I,2)'-fc)/(1-fc)]';%normalizing
    end
  else
    try
      %if isa(G,'fdata') || isa(G,'struct')
      Gi = G;
      cdf = ['cdf',getdistname(Gi.distribution)];
      fc = feval(cdf,c,Gi);
      G = [z(:),(feval(cdf,z(:),Gi)-fc)/(1-fc)];
    catch
      error('Wrong input for G! Must be a 2 column matrix')
    end
  end
  if options.wdata
    G = createwdata('data',G(:,2:end),'args',G(:,1),'labels',{xlab,ylab},'note','EDFCND');  
    
    G = wdata(G);
    
  end
end

if options.wdata
  Fz = createwdata('data',Fz1,'args',z,'labels',{xlab,ylab},'note','EDFCND');
  Fz = wdata(Fz);
else
  Fz = [z(:),Fz1(:)];
end
  F = Fz;

