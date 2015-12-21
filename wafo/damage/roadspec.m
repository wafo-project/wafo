function S = roadspec(sdata,a,C)
%ROADSPEC Spectral density (frequency) for a road 
%
% CALL:  S = roadspec(data,a,C);
%
%  Output:
%        S = the spectral density (structure array)
%  Input: 
%    sdata = the data vector [wl wu n], where 
% 
%	wl = lower truncation frequency  (default 4/257)
%	wu = upper truncation frequency  (default 4)
%	n  = number of evaluation points (default 257)
%     a,C  = constants in the spectral density
%
% The model is given by 
% 
%  	              S(w) = C/(w^a),  wl < w < wu
% 
% Usually 2 < a < 3, see the literature. For the value of  c, 
% Kamash and Robson (1978) give the values
%
% Motorway:  3e-8 < C < 50e-8
% Highway:   3e-8 < C < 800e-8
% Minor highway: 50e-8 < C < 3000e-8
% 

% References:
% Lindgren, G. (1981).
% Jumps and Bumps on Random Roads.
% Journal of Sound and Vibration, Vol 78, pp 383-395
%
% Kamash, K.M.A., and Robson, J.D. (1978).
% The Application of Isotropy in Road Surface Modelling 
% Journal of Sound and Vibration, Vol 57, pp 89-100
%
% Jogréus, C. (1983).
% Fordonsrörelser och stokastiska vägmodeller.
% Master's thesis, Mathematical Statistics, Lund University 

% Tested on Matlab 6.0
% History:
% Modified by jr 01-April-2001
% - structure array introduced
% - help text modified
% By Mats Frendahl, 1993


if nargin<1 || isempty(sdata), sdata = [4/257 4 357]; end 
if nargin<2 || isempty(a), a = 2.1; end 
if nargin<3 || isempty(C), C = 50e-8; end
if (a<2) || (a>3)
  disp(' The parameter  a  must be in (2,3). Program will terminate.')
  return
end

wl = sdata(1); wu = sdata(2); n = sdata(3); 
wv = linspace(0,wu,n);
% wv = linspace(w0:(l1-l0)/(n-1):l1;  % OLD freq vector
spv = C*(wv.^a).^(-1);

S=createspec;
S.S=spv; 
S.w=wv;
S.type='freq';
S.note='Spectrum: road spectrum';
S.S(wv<wl)=0;
S.S(1)=0; % must be zero at zero freq since discrete spectrum

