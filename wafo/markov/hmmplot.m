function hmmplot(varargin)
% HMMPLOT  plots a Hidden Markov Model.
%
% CALL: hmmplot(x,z)
%       hmmplot(x,z,t,z_range,Title,Ylabel,colType,fontsize)
%
%   x       = Switching process.
%   z       = Regime process.
%   t       = Time.
%   z_range = Range of regime, e.g. [1 3] for 3 regime states.
%   Title   = Title of plot.
%   Ylabel  = Y-label of x-plot.
%   colType = 0: One colour (default), 1: Different colours for each regime
%   fontsize= Fontsize of text.
%
% Examples: Switching AR(1)-process. (Example 3 in thesis)
%   P= [0.975 0.02 0.005; 0.01 0.98 0.01; 0.005 0.02 0.975];
%   C = [1 1 1]';
%   A = [1 -0.5; 1 -0.3; 1 0.5];
%   m = [-1 0 3]';
%   s2 = [1 1 1.44]';
%   [x,z] = sarmasim(C,A,m,s2,P,500);
%   hmmplot(x,z)
%   hmmplot(x,z,(0:499)/400,[1 3],'Switching AR(1)-process','X(t)')
%   hmmplot(x,z,[],[1 3],'','',1,20)

% Copyright (c) 1997 by Pär Johannesson
% Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

plothmm(varargin{:})
end
