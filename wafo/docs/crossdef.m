% CROSSDEF level v crossing definitions and nomenclature 
%
% Definition of level v crossing:   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Let the letters 'm', 'M', 'F', 'L','d' and 'u' in the 
% figure below denote local minimum, maximum, first value, last 
% value, down- and up-crossing, respectively. The remaining 
% sampled values are indicated with a '.'. Values that are identical 
% with v, but do not cross the level is indicated with the letter 'o'.
%
% We have a level up-crossing at index, k, if 
%
%           x(k) <  v and v < x(k+1) 
% or if 
%           x(k) == v and v < x(k+1) and x(r) < v for some di < r <= k-1 
%
% where di is  the index to the previous downcrossing.
% Similarly there is a level down-crossing at index, k, if 
% 
%           x(k) >  v and v > x(k+1)
%  or if 
%           x(k) == v and v > x(k+1) and x(r) > v  for some ui < r <= k-1 
%
% where ui is  the index to the previous upcrossing.
%
% The first (F) value is a up crossing if  x(1) = v and x(2) > v. 
% Similarly, it is a down crossing if      x(1) = v and x(2) < v. 
%
%      M
%    .   .                  M                   M
%  .      . .             .                   .   . 
%F            d               .             .       L
% ----------------------u-------d-------o--------------------- level v
%               .     .           .   .   u
%                 .                 m
%		   m
%
% See also  perioddef, wavedef, tpdef, findcross, dat2tp

% history:
% by Per A. Brodtkorb 19.09.1999

more on
help crossdef
more off
