function ltype = lagtype(R)
%LAGTYPE Returns the lag type of a Covariance struct.
%
% CALL:  ltype = lagtype(R)
%
%  ltype = Character vector containing:
%          'x' if lag of first space dimension is given.
%          'y' if lag of second space dimension is given.
%          't' if time lag is given
%      R = Covariance function structure
%
% Example:
%  R = spec2cov(jonswap);
%  lagtype(R)
%
% See also  datastructures

% Tested on: matlab 5.2
% History:
% by pab 11.10.2001

names=fieldnames(R); 
ind=find(strcmp(names,'x')+strcmp(names,'y')+strcmp(names,'t')); %options are 'x' and 'y' and 't' 
if isempty(ind)
  error('This is not a Covariance structure')
 % ltype=[];
else
  ltype=char(names(ind)).'; 
end
