function [varargout]= snplot(varargin)
%SNPLOT Plots SN-data and estimates parameters 
% according to the least-squares method
%
% CALL:  [e,beta,s2] = snplot(S,N,nr);
%
%  where
%
%        e,beta = the parameters in the model, 
%        s2     = the residual variance for linear regression 
%                 in logS/logN plane,
%        S      = a nx1 vector with S-data,
%        N      = a nx1 vector with N-data,
%        nr     = plot parameter (optional input argument);
%
%                 1 = only SN-data will be plotted,
%                 2 = SN-data and fitted line will be plotted,
%                 3 = only log(S)/log(N)-data will be plotted,
%                 4 = log(S)/log(N)-data and fitted line will be plotted,
%                11-14 = Same as above but x-axis = N, y-axis = S . 
%  
% Model:
%       N(s) = K/(e*s^beta)
%
% Example:  
%   sn = load('sn.dat'); s = sn(:,1); N = sn(:,2);
%   [e,beta,s2] = snplot(s,N,2);   % S-N, x-axis = S, y-axis = N
%   [e,beta,s2] = snplot(s,N,12);  % N-S, x-axis = N, y-axis = S
%
%   close all;
  
% Tested on: Matlab 6.0
% History:
% Correction by PJ 07-Jul-2005
%   Changed 'break' to 'return'
% Revised by jr 01-Apr-2001
% - example, help text
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Original version from FAT by Mats Frendahl 
%   Copyright 1993, Mats Frendahl, Dept. of Math. Stat., University of Lund.

[varargout{1:nargout}] = plotsn(varargin{:});
return