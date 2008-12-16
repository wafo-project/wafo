function [x,out2,out3] = mctpsim(varargin) 
%MCTPSIM  Simulates a Markov chain of turning points
%  The process x has the state space {1,2,...,n}.
%
% CALL: [x] = mctpsim(F,T);
%       [x] = mctpsim(F,T,init);
%       [x] = mctpsim(F,T,init,'x');
%       [RFM,RFM0,res] = mctpsim(F,T,init,'RFM');
%       [x,RFM] = mctpsim(F,T,init,'x,RFM');
%
% x       = Simulated switching Markov turning points.
% RFM     = Rainflow matrix for x.                        [nxn]
% RFM0    = Rainflow matrix for x (without the residual). [nxn]
% res     = Residual from rainflow count.                 [nx2]
%
% F       = Cell array of min-Max and Max-min matrices {1,2}
% F{1,1}  = min-Max matrix, process 1                  [nxn]
% F{1,2}  = Max-min matrix, process 1                  [nxn]
% T       = Length of simulation.
% init.x0 = Initial state of process x. If not given, it will start from
%           the stationary distribution of minima.
%
% Simulates a Markov chain of turning points,
% The process x has the state space {1,2,...,n}.
%
% If a matrix F{1,2}=[], then the process will
% be assumed to be time-reversible.
%
% Examples: 
%   FF = mktestmat([-1 1 32],[-0.2 0.2],0.15,1);
%   x = mctpsim({FF []},1000); 

init.z0=1;

[x,out2,out3] = smctpsim(1,varargin{:});
%switch whatOut

%case {'x'}
  
%  x = smctpsim(1,vargin{:});

%case {'RFM'}
  
%  [x,out2,out3] = smctpsim(1,vargin{:});

%case {'x','x,RFM'}

%  [x,out2] = smctpsim(1,vargin{:});

%end
