function figsort(varargin)
%FIGSORT Sort figure windows in ascending order.
%
%   FIGSORT brings up the figures in ascending order so the figure with the
%   lowest figure number becomes the first window.
%   FIGSORT(FIGS) can be used to specify which figures that should be
%   sorted.
%
% See also  figcycle, figkeep, figmaximize, figrestore, figpile, figstack,
%          figtile

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-06-07 13:30:52
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://home.online.no/~pjacklam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the handles to the figures to process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figs = [];
for ix=1:nargin
  currArg = varargin{ix};
  if isnumeric(currArg) 
    figs = [figs,currArg];
  else
    tmp = double(currArg);
    try
      if any( double('0')<= tmp & tmp<= double('9') ) || strcmpi(currArg,'gcf')
        figs = [figs, eval(currArg)];
      end
    end
  end
end
if isempty(figs)                          % If no input arguments...
   figs = findobj('Type', 'figure');    % ...find all figures.
   figs  = sort(figs);
end

if isempty(figs)
   disp('No open figures or no figures specified.');
   return
end

nfigs = length(figs);

for i = nfigs:-1:1
   figure(figs(i));             % Bring up i'th figure.
end
