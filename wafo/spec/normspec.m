function [varargout]=normspec(varargin)
%NORMSPEC Normalize a spectral density such that m0=m2=1
%
% CALL:  [Sn,mn4,m0,m2,m4,m1]=normspec(S,m0,m2,plotflag);
%
%        Sn  = the normalized spectrum struct
%        mn4 = the fourth spectral moment of Sn.
%        mi  = the i'th spectral moment of S.
%        S   = the original spectrum struct
%   plotflag = 0, do not plot the normalized spectrum (default).
%              1, plot the normalized spectrum.
%           
% Normalization performed such that
%    INT S(freq) dfreq = 1       INT freq^2  S(freq) dfreq = 1
% where integration limits are given by  freq  and  S(freq)  is the 
% spectral density; freq can be frequency or wave number.
% The normalization is defined by
% A=sqrt(m0/m2); B=1/A/m0; freq'=freq*A; S(freq')=S(freq)*B;
%
% If S is a directional spectrum then a normalized gravity (.g) is added
% to Sn, such that mxx normalizes to 1, as well as m0 and mtt.
% (See spec2mom for notation of moments)
%
% If S is complex-valued cross spectral density which has to be
% normalized, then m0, m2 (suitable spectral moments) should be given.
%
% Example: 
%   S = jonswap;
%   [Sn,mn4] = normspec(S);
%   assert(spec2mom(Sn,2), [1,1], 1e-4)     % Should be equal to one!
% 
% See also  spec2mom, spec2spec 

% Tested on: Matlab 5.3
%
% History:
% revised ir 31.08.01, introducing normalizing  g for w-spectrum
% Revised by jr 20.04.2001
% - Condition on nargout related to calculation of mn4; minor change
% - Updated information, added example
% By es 28.09.1999


[varargout{1:nargout}]=specnorm(varargin{:});

