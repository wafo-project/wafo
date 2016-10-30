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
%  assert(S.S(100:105), [2.18098263255750, 2.38730066654042,...
%                           2.59856918262022, 2.81441896250507,...
%                           3.03482079647815, 3.26016190326422]', 1e-9)
%  assert(S.w(100:105), [ 0.440859375000000, 0.445312500000000,...
%                         0.449765625000000, 0.454218750000000,...
%                         0.458671875000000, 0.463125000000000]', 1e-9)
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

