% TPDEF turning points definitions and numenclature
%
% Definition of turningpoints: 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Local minimum or maximum are indicated with the
% letters 'm' or 'M'. Turning points in this connection are all
% local max (M) and min (m) and the last (L) value and the 
% first (F) value if the first local extremum is a max. 
%
% (This choice is made in order to get the exact up-crossing intensity
% from rfc by mm2lc(tp2mm(rfc)) )
%
% Definition of min2min and Max2Max wave:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% A min2min wave (mw) is defined starting from a minimum and  
% ending in the following minimum.
% Similarly a Max2Max wave (Mw) is thus a wave from a maximum 
% to the next maximum (all waves optionally rainflow filtered).
%
% Definition of Max2min and min2Max cycle:
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% A min2Max cycle (mM) is defined as the pair of a minimum
% and the following Maximum. Similarly a Max2min cycle (Mm) 
% is defined as the pair of a Maximum and the following minimum.
% (all cycles possibly rainflowfiltered)
%
%              <----- Direction of wave propagation 
%
%   M                  M
%  / \       .... M   /:\_           M_            M
% F   \     |    / \m/ :  \         /: \          / \  
%      \  h |   /      :   \       / :  \        /   \
%       \   |  /       :    \     /  :   \_    _/     \_   L
%        \_ | /        :     \_m_/   :     \m_/         \m/
%          \m/         :             :      :            :
%                      <------Mw----->      <-----mw----->
%
% See also  wavedef, crossdef, dat2tp

% history:
% by Per A. Brodtkorb 19.09.1999
% revised pab 15.05.2000


more on
help tpdef
more off
