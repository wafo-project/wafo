function [varargout] = welch_cpsd(varargin)
%WELCH_CPSD Estimate cross power spectral density using Welch's method
%
% CALL: [Pxy, f, options, Pci] = welch_cpsd(x, y, options)
%
% WELCH_CPSD estimate cross spectrum density of a pair of signals. This chops the 
% signals into overlapping sections, windows each section and applies a Fourier
% transform to determine the frequency components at that slice. The
% magnitudes of these sections are then averaged to produce the estimate Pxy.
% The confidence interval around the estimate is returned in Pci.
%
% See welch_psd for an explanation of the available parameters.
%
% Note welch_cpsd(x,y) == conj(csd(x,y))*2/Fs = csd(y,x)*2/Fs
%
% Example
%  Fs = 400; x = 0:1/Fs:1; w1 = 2*pi*30; w2 = 2*w1;
%  y1 = cos(w1*x) + randn(size(x));
%  y2 = cos(w2*x) + randn(size(x));
%  [Pxy,fi] = welch_cpsd(y1,y2);
%  plot(fi,Pxy)
%
% See also welch_psd, welch_tfe, welch_cohere


% History
% Revised pab april 2007
% - renamed from csd to welch_cpsd
% - updated help header to conform to the wafo style.

% Copyright (C) 2000 Paul Kienzle
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  if nargin < 2
    disp('Pxy=welch_cpsd(x,y,...)  [see welch_psd for details]'); 
    return
  end
  if nargout==0, 
    welch_psd(varargin(1:2),'cpsd',varargin{3:nargin});
  else
    [varargout{1:nargout}] = welch_psd(varargin(1:2),'cpsd',varargin{3:nargin});
  end
end
