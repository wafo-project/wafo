%WAVEDEF Wave definitions and nomenclature 
%
% Definition of trough and crest: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% A trough (t) is defined as the global minimum between a
% level v down-crossing (d) and the next up-crossing (u)  
% and a crest (c) is defined as the global maximum between a
% level v up-crossing and the following down-crossing.
%
% Definition of down- and up -crossing waves: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% A level v-down-crossing wave (dw) is a wave from a 
% down-crossing to the following down-crossing.
% Similarly, a level v-up-crossing wave (uw) is a wave from an up-crossing 
% to the next up-crossing.
%
% Definition of trough and crest waves: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% A trough-to-trough wave (tw) is a wave from a trough (t) to the 
% following trough. The crest-to-crest wave (cw) is defined similarly.
% 
% Definition of min2min and Max2Max wave:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% A min2min wave (mw) is defined starting from a minimum (m) and  
% ending in the following minimum.
% Similarly a Max2Max wave (Mw) is thus a wave from a maximum (M) 
% to the next maximum (all waves optionally rainflow filtered).
%
%           <----- Direction of wave propagation
%
%
%   <------Mw-----> <----mw---->
%   M             : :  c       :
%  / \            M : / \_     :     c_            c
% F   \          / \m/    \    :    /: \          /:\  
%------d--------u----------d-------u----d--------u---d-------- level v
%       \      /:           \  :  /: :  :\_    _/  : :\_   L
%        \_   / :            \_t_/ : :  :  \t_/    : :  \m/
%          \t/  <-------uw---------> :  <-----dw----->
%           :                  :     :             : 
%           <--------tw-------->     <------cw----->              
%
% (F= first value and L=last value).
%
% See also  tpdef, crossdef, dat2tc, dat2wa, dat2crossind

% history:
% revised pab 21.01.2000
% -  changed  from: Wave propagating direction 
%               to: Direction of wave propagation
% by Per A. Brodtkorb 19.09.1999

more on
help wavedef
more off
