function S=demospec(stype)
% DEMOSPEC Loads a precreated spectrum of chosen type
%
% CALL:  S=demospec(stype)
%
%       S = spectrum struct
%   stype = string 'freq' or 'dir' (default 'freq')
%
% Spectra of types 'freq' and 'dir' are available
% and can be used as test or demo examples.
% Any other type is possible to get using SPEC2SPEC
% 
% Example
%  S = demospec('freq');
%  plotspec(S)
%
% See also  createspec, datastructures, spec2spec

% Tested on: Matlab 5.3
% History: 
% revised by es 05.06.00 changed call to SPREADING since spreading changed
%  by es 21.12.1999 m-file calls instead of mat-files
%  by es 23.09.1999

if nargin<1||isempty(stype)
  stype='freq';
end

if strcmpi(stype,'freq')
  S=jonswap(linspace(0,1.14,257));
else
  S=jonswap(linspace(0,1.14,257));
  D=spreading(linspace(-pi,pi,101),'cos2s',0,15,S.w,0);
  S=mkdspec(S,D);
end
S.note=['Demospec: ',S.note,', truncated at 2*wp'];

