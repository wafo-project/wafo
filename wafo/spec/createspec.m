function S=createspec(stype,freqtype)
% CREATESPEC Spectrum struct constructor
%
% CALL:  S=createspec(stype,freqtype)
% 
%      stype = 'freq'   Frequency spectrum (default)
%              'dir'    Directional spectrum
%              'k1D'    Wave number spectrum 1D
%              'k2D'    Wave number spectrum 2D
%              'encdir' Encounter directional spectrum
%              'enc'    Encounter frequency spectrum
%
%   freqtype = 'w' angular frequency (rad/sec) (default)
%              'f' frequency         (Hz)
%
% Example: Create a struct with proper fieldnames for directional spectrum
%    S=createspec('dir');
%    S.date = '';
%    names = {'S', 'w', 'theta', 'tr', 'h', 'type', 'phi','norm','note', 'date'};
%    vals = {[],[],[],[],inf,'dir', 0, 0, [], ''};
%    assert(fieldnames(S), names');
%    assert(struct2cell(S), vals');
% 
% See also  createcov, datastructures

% Tested on: Matlab 5.3
% History:
%  revised by IR 03.04.2001 - added S.phi=0.
%  revised by jr 10.07.2000 - Line 47 and 49: semicolon added
%  revised by es 25.05.2000 - help-text changes  
%  revised by jr 14.01.2000 - Field added: norm: 0 
%  revised by es 19.09.1999
%  by pab 12.08.99

if nargin<1||isempty(stype)
  stype='freq';
else
  stype=lower(stype);
end

if nargin<2||isempty(freqtype)
  freqtype='w';
else
  freqtype=lower(freqtype);
end

n=length(stype);
S=struct('S',[]);
if strcmp(stype(1),'k') % wavenumber spectrum
  S.k=[];
  if strcmp(stype(n-1),'2')
    S.k2=[];
  end
else % 
  if strcmp(freqtype,'f')
    S.f=[];
  else
    S.w=[];
  end
  if strcmp(stype(n-2:n),'dir')
   S.theta=[];
  end
end

S.tr=[];

if strcmp(stype(1:3),'enc')
   S.v=0;
   S.phi=0;
elseif strcmp(stype(1:3),'rot')
   S.phi=0;
end
S.h=inf;
S.type=stype;
S.phi=0.;
S.norm=0;
S.note=[];
S.date=datestr(now);
