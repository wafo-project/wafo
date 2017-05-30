function figkeep(varargin)
%FIGKEEP Keeps figure windows of your choice and clears the rest.
%
%  CALL figkeep(figs)
%  
%  figs = vector of integers specifying which figures to keep. 
% 
% Example: %keep only figures 1,2,3,5 and 7
%  for ix = 1:10,figure(ix),end
%  figkeep 1:3  5 7 
%  % or 
%  figkeep([1:3  5 7])  
%
%  close all;
%  
% See also close, figcycle, figmaximize, figrestore, figpile, figstack,
%          figsort, figtile
  
% by Per A. Brodtkorb 05 Apr 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the handles to the figures to process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figs2keep = [];
for ix=1:nargin
  currArg = varargin{ix};
  if isnumeric(currArg) 
    figs2keep = [figs2keep,currArg];
  else
    tmp = double(currArg);
    try
      if any( double('0')<= tmp & tmp<= double('9') ) || strcmpi(currArg,'gcf')
        figs2keep = [figs2keep, eval(currArg)];
      end
    end
  end
end

if isempty(figs2keep)
   disp('No figures specified to keep.');
   return
end
allfigs = findobj('Type', 'figure');    % ...find all figures.
% Remove figure handles in the "keep" list
figs2delete = setdiff(allfigs,figs2keep);
if size(figs2delete,1)>1,
  figs2delete = figs2delete.';
end
close(figs2delete)
return