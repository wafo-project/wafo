% PERIODDEF wave periods (lengths) definitions and nomenclature
%
% Definition of wave periods (lengths):
%---------------------------------------
%
%
%           <----- Direction of wave propagation 
%
%               <-------Tu--------->
%               :                  :
%               <---Tc----->       :
%               :          :       : <------Tcc---->
%   M           :      c   :       : :             :
%  / \          : M   / \_ :       : c_            c
% F   \         :/ \m/    \:       :/  \          / \  
%------d--------u----------d-------u----d--------u---d-------- level v
%       \      /            \     /     :\_    _/:   :\_   L
%        \_   /              \_t_/      :  \t_/  :   :  \m/
%          \t/                 :        :        :   :
%           :                  :        <---Tt--->   : 
%           <--------Ttt------->        :            :
%                                       <-----Td----->
%  Tu   = Up crossing period
%  Td   = Down crossing period
%  Tc   = Crest period, i.e., period between up crossing and 
%         the next down crossing 
%  Tt   = Trough period, i.e., period between down crossing and 
%         the next up crossing
%  Ttt  = Trough2trough period
%  Tcc  = Crest2crest period
%
%
%           <----- Direction of wave propagation 
%
%                <--Tcf->                              Tuc
%                :      :               <-Tcb->        <-> 
%   M            :      c               :     :        : :
%  / \           : M   / \_             c_    :        : c
% F   \          :/ \m/    \           /  \___:        :/ \  
%------d---------u----------d---------u-------d--------u---d-------- level v
%      :\_      /            \     __/:        \_    _/     \_   L
%      :  \_   /              \_t_/   :          \t_/         \m/
%      :    \t/                 :     :
%      :     :                  :     : 
%      <-Ttf->                  <-Ttb-> 
%                
% 
%  Tcf  = Crest front period, i.e., period between up crossing and crest
%  Tcb  = Crest back period, i.e., period between crest and down crossing 
%  Ttf  = Trough front period, i.e., period between down crossing and trough
%  Ttb  = Trough back period, i.e., period between trough and up crossing 
%
% Also note that Tcf and Ttf can also be abbreviated by their crossing
% marker, e.g. Tuc (u2c)  and Tdt (d2t), respectively. Similar applies to all the
% other wave periods and wave lengths. 
% (The nomenclature for wave length is similar, just substitute T and
% period with L and length, respectively)
%
%             <----- Direction of wave propagation 
%
%                      <--TMm-->
%           <-TmM->    :       :
%   M       :     :    M       :
%  / \      :     M   /:\_     :     M_            M
% F   \     :    / \m/ :  \    :    /: \          / \  
%      \    :   /      :   \   :   / :  \        /   \
%       \   :  /       :    \  :  /  :   \_    _/     \_   L
%        \_ : /        :     \_m_/   :     \m_/         \m/
%          \m/         :             :      :            :
%                      <-----TMM----->      <----Tmm----->
%
%
%  TmM = Period between minimum and the following Maximum
%  TMm = Period between Maximum and the following minimum
%  TMM = Period between Maximum and the following Maximum
%  Tmm = Period between minimum and the following minimum
%
% See also  wavedef, ampdef, crossdef, tpdef

% history:
% revised pab 15.05.2000
% - minor corrections for Tcb
% - added Tmm, TMM, etc.....
% revised pab 21.01.2000
% - updated help header
% -  changed from: Wave propagating direction 
%              to: Direction of wave propagation
% by Per A. Brodtkorb 19.09.1999


more on
help perioddef
more off
