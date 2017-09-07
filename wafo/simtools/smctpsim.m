function [x,z,TT] = smctpsim(P,F,T,init,whatOut)
%SMCTPSIM  Simulates a switching Markov chain of turning points,
%   i.e. a switching process with a Markov chain of turning points
%   within each regime.
%   The switching process x has the state space {1,2,...,n} and
%   the regime process z has the state space {1,2,...,r}.
%
% [x,z] = smctpsim(P,F,T);
% [x,z] = smctpsim(P,F,T,init);
% [x,z] = smctpsim(P,F,T,init,'x');
% [RFM,RFM0,res] = smctpsim(P,F,T,init,'RFM');
% [x,z,RFM] = smctpsim(P,F,T,init,'x,RFM');
%
% x       = Simulated switching Markov turning points.
% z       = Simulated regime process.
% RFM     = Rainflow matrix for x.                        [n,n]
% RFM0    = Rainflow matrix for x (without the residual). [n,n]
% res     = Residual from rainflow count.                 [n,2]
%
% P       = Transition matrix for regime process.      [r,r]
% F       = Cell array of min-Max and Max-min matrices {r,2}
% F{i,1}  = min-Max matrix, process i                  [n,n]
% F{i,2}  = Max-min matrix, process i                  [n,n]
% T       = Length of simulation.
% init.x0 = Initial state of process x. If not given, it will start from
%          the stationary distribution of minima given z(1).
% init.z0 = Initial state of regime process. If not given, it will start 
%          from the stationary distribution of the Markov chain.
%
% If a matrix F{i,2}=[], then the process will
% be assumed to be time-reversible.
%
% One can also simulate a switching process other than Markov switching.
% In that case P is a struct array describing the switching
% P.P        = Embedded transition matrix
% P.distr    = 'exp':       Exponential distribution (Markov switching)
%              'phasetype': Phasetype distribution
%              'const':     Constant
%              'other':     Other distribution (simTfun)
% P.simTfun  = Function for simulating occupation times.
% P.Tpar     = Parameters to distributions       {r,1}
% P.sequence = Specified sequence.               [qx2]
%
% Examples: Two regime states.
%  [x,z] = smctpsim(P,F,T); % Two regime states
% Example: A Markov chain of turning points (One regime state).
%  [x,z] = smctpsim(1,F,T);

% Tested  on Matlab  5.3
%
% History:
% Correction by PJ 20-Jan-2004
%   Some corrections concerning initial state of regime process
%   and simulation of non-Markovian regime process.
% Correction by PJ 12-Aug-2002
%   'RFM' opt didn't work cause x was not initiated. Now fixed!
% Revised by PJ 19-May-2000
%   Corrected method for simulating starting conditions.
% Revised by PJ Jan-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1997
%   Copyright (c) 1997 by Pär Johannesson
%   Toolbox: Rainflow Cycles for Switching Processes V.1.0, 2-Oct-1997

TT(1,:) = clock;

ni = nargin;
%no = nargout;
%error(nargchk(3,5,ni));
narginchk(3,5)
Zstr = '123456789';

if ni < 4,  init = []; end
if ni < 5,  whatOut = []; end

if isempty(init)
    init.x0 = [];
    init.z0 = [];
end
if ~isfield(init,'x0'), init.x0=[]; end
if ~isfield(init,'z0'), init.z0=[]; end

if isempty(whatOut)
    whatOut = 'x';
end

if isa(P,'double')
    Ptype = 'P';
elseif isa(P,'struct')
    Ptype = 'struct';
    S = P;
else
    error('P should be matrix or struct-array.');
end

r = size(F,1);   % Number of regime states
%r = length(P);   % Number of regime states

n = length(F{1,1});  % Number of levels

% Check that the rowsums of P are equal to 1

if strcmp(Ptype,'struct')
    if isfield(S,'P')
        P=S.P;
    else
        P=[];
    end
end

if ~isempty(P)
    sumP = sum(P,2).';
    if sum(sumP == 1) ~= length(P)
        warning('WAFO:SMCTPSIM',' Rowsums of P not equal to 1. Renormalizing!');
        for i = 1:length(P)
            P(i,:) = P(i,:)/sumP(i);
        end
    end
end

TT(2,:) = clock;

% Normalize the rowsums of F{1,1},...,F{r,1} to 1
%  ==>  Q{1,1},...,Q{r,1}

for i = 1:r
    QQ{i,1} = triu(F{i,1},1); % Set zeros below diagonal and diagonal
    % Set negative elements to zero
    [I,J] = find(QQ{i,1}<0);
    if length(I) ~= 0
        warning('WAFO:SMCTPSIM',['Negative elements in Q' Zstr(i) '. Setting to zero!']);
        for k = 1:length(I)
            QQ{i,1}(I(k),J(k)) = 0;
        end
    end
    
    sumQQi = sum(QQ{i,1},2).';
    % Normalize rowsums
    if sum(sumQQi == 1) ~= length(QQ{i,1})
        %disp(['Warning: Rowsums of Q' Zstr(i) ' not equal to 1. Renormalizing!']);
        for j = 1:n-1
            if sumQQi(j)~=0, QQ{i,1}(j,:) = QQ{i,1}(j,:)/sumQQi(j); end
        end
    end
end

TT(3,:) = clock;

% Normalize the rowsums of F{1,2},...,F{r,2} to 1
%  ==>  Q{1,2},...,Q{r,2}

% Normalize the rowsums of Fh1,...,Fhr to 1  ==>  Q1,...,Qr

for i = 1:r
    if isempty(F{i,2})        % Time-reversible
        QQ{i,2} = F{i,1}';
    else
        QQ{i,2} = F{i,2};
    end
    
    QQ{i,2} = tril(QQ{i,2},-1); % Set zeros above diagonal and diagonal
    % Set negative elements to zero
    [I,J] = find(QQ{i,2}<0);
    if length(I) ~= 0
        warning('WAFO:SMCTPSIM',['Negative elements in Qh' Zstr(i) '. Setting to zero!']);
        for k = 1:length(I)
            QQ{i,2}(I(k),J(k)) = 0;
        end
    end
    
    sumQQi = sum(QQ{i,2},2).';
    if sum(sumQQi == 1) ~= length(QQ{i,2})
        %disp(['Warning: Rowsums of Qh' Zstr(i) ' not equal to 1. Renormalizing!']);
        for j = 2:n
            if sumQQi(j)~=0, QQ{i,2}(j,:) = QQ{i,2}(j,:)/sumQQi(j); end
        end
    end
end

TT(4,:) = clock;



% Initial state of z0, for regime process
% and x0, for X-process, start from a minimum

% Make the transition matrix Q for the joint MC (X_k,Z_k)
[Q] = smctp2joint(P,QQ);

% Stationary distribution = ro of minima
[ro_min,ro_max] = mctp2stat(Q);  
ro = ro_min;

% Start values
e0 = rand(1,1);  % Generate random numbers
if isempty(init.z0) && isempty(init.x0)
    x0z0 = min(find( e0<=cumsum(ro) ));
    x0 = floor((x0z0+1)/r);
    z0 = mod(x0z0-1,r)+1;
elseif isempty(init.x0)
    z0 = init.z0;
    rox0 = ro(z0:r:end); % Pick stat. distr. for regime z0
    rox0 = rox0/sum(rox0);
    x0 = min(find( e0<=cumsum(rox0) ));
elseif isempty(init.z0)
    x0 = init.x0;
    z0 = [];  % Start from stat. distr of P
else % Both z0 znd x0 are given
    x0 = init.x0;
    z0 = init.z0;
end


if strcmp(Ptype,'struct')
    if isfield(S,'P')
        S.P=P;
    end
end


%
% Initiate vectors
%

switch whatOut
    
    case {'x','x,RFM'}
        if strcmp(Ptype,'P')
            z = mcsim(P,T,z0);  % Simulate Regime
        else
            z=zeros(T,1);
            z(1)=init_z(S,init);
        end
        e=rand(T,1);
        x=zeros(T,1);
        
    case {'RFM'}
        if strcmp(Ptype,'P')
            z = mcsim(P,1,z0);  % Simulate Regime
        else
            z=init_z(S,init);
        end
        e=rand(1,1);
        
end

x(1) = x0;

TT(5,:) = clock;

%
% Simulate Switching Markov turning points
%

% Calculate cumulative distributions

cumsumP = cumsum(P,2);
ones_n = ones(n,1);
for i = 1:r 
    cumsumQQ{i,1} = cumsum(QQ{i,1},2);
    cumsumQQ{i,2} = cumsum(QQ{i,2},2);
    cumsumQQ{i,1}(:,end) = ones_n;
    cumsumQQ{i,2}(:,end) = ones_n;
end

% Simulate

switch whatOut
    
    case {'x','x,RFM'}
        switch Ptype
            case 'P'
                for k=2:T
                    if rem(k,2) == 0 % min-to-Max
                        x(k) = sum( cumsumQQ{z(k),1}(x(k-1),:) <= e(k) ) + 1;
                    else             % Max-to-min
                        x(k) = sum( cumsumQQ{z(k),2}(x(k-1),:) <= e(k) ) + 1;
                    end
                end
                
            case 'struct'
                t0 = simT(S,z(1));
                for k=2:T
                    if t0 > 1
                        z(k) = z(k-1);
                        t0=t0-1;
                    else
                        z(k-1:k) = mcsim(S.P,2,z(k-1));
                        t0 = simT(S,z(k));
                    end
                    %      fprintf(1,'k=%d, z=%d\n',k,z(k));
                    if rem(k,2) == 0 % min-to-Max
                        x(k) = sum( cumsumQQ{z(k),1}(x(k-1),:) <= e(k) ) + 1;
                    else             % Max-to-min
                        x(k) = sum( cumsumQQ{z(k),2}(x(k-1),:) <= e(k) ) + 1;
                    end
                end
        end
        
        
    case 'test' %'x,RFM'
        [RFM0,res,nres] = dtp2rfm_init(n);
        [RFM0,res,nres] = dtp2rfm1(x(1),RFM0,res,nres);
        for k=2:T
            e=rand(2,1);
            
            % Simulate Regime
            z(k) = sum( cumsumP(z(k-1),:) <= e(1) ) + 1;
            
            % Simulate MCTP
            if rem(k,2) == 0 % min-to-Max
                x(k) = sum( cumsumQQ{z(k),1}(x(k-1),:) <= e(2) ) + 1;
            else             % Max-to-min
                x(k) = sum( cumsumQQ{z(k),2}(x(k-1),:) <= e(2) ) + 1;
            end
            [RFM0,res,nres] = dtp2rfm1(x(k),RFM0,res,nres);
        end
        [RFM] = dtp2rfm_collect(RFM0,res,nres);
        
    case 'RFM'
        [RFM0,res,nres] = dtp2rfm_init(n);
        [RFM0,res,nres] = dtp2rfm1(x,RFM0,res,nres);
        for k=2:T
            e=rand(2,1);
            
            % Simulate Regime
            z = sum( cumsumP(z,:) <= e(1) ) + 1;
            
            % Simulate MCTP
            if rem(k,2) == 0 % min-to-Max
                x = sum( cumsumQQ{z,1}(x,:) <= e(2) ) + 1;
            else             % Max-to-min
                x = sum( cumsumQQ{z,2}(x,:) <= e(2) ) + 1;
            end
            [RFM0,res,nres] = dtp2rfm1(x,RFM0,res,nres);
        end
        [RFM] = dtp2rfm_collect(RFM0,res,nres);
        
end

TT(6,:) = clock;

% Output arguments

switch whatOut
    case 'x,RFM'
        RFM = dtp2rfm(x,n);
        TT = RFM;
    case 'RFM'
        x = RFM;
end


function z=init_z(S,init)

% Initiate regime.

if isempty(init.z0)
  r=length(S.P);
  m = zeros(1,r);
  switch S.distr
    case {'exp','const'},
      for i=1:r, m(i)=S.Tpar{i}; end
    case 'phasetype',
      for i=1:r, m(i)=sum(S.Tpar{i}); end
  end
  statP=mc2stat(S.P);
  ro = statP.*m;
  z = sum( cumsum(ro) <= rand ) + 1;
else
  z=init.z0;
end

    
    function t=simT(S,z)
    
    % Simulate occupation time.
    
    if iscell(S.distr)
        distr = S.distr{z};
    else
        distr = S.distr;
    end
    
    switch distr
        
        case 'exp'
            t = simFS(S.Tpar{z});
        case 'phasetype'
            t = simPhaseType(S.Tpar{z});
        case 'const'
            t = S.Tpar{z};
        case 'other'
            error('other: Not yet implemented.');
    end
    
    % Simulate First Success R.V.
    
    function t=simFS(m)
    % Simulation method taken from stats-toolbox 'geornd'
    
    u = rand(size(m));
    p = 1./m;
    t = floor(log(u) ./ log(1 - p)) + 1;
    %t = - m .* log(rand(size(m)));
    
    % Simulate Exponential R.V.
    function t=simExp(m)
    
    t = - m .* log(rand(size(m)));
    
    % Simulate discrete Phasetype R.V.
    
    function t=simPhaseType(m)
    
    t = sum(simFS(m));
    %t = sum(simExp(m));
    
    % Description of variables for  dtp2rfm
    %
    % RFM   = Rainflow Matrix (residual included).    [nxn]
    % RFM0  = Rainflow matrix (without resudual).     [nxn]
    % res   = Residual.                               [2*n,1]
    % nres  = Length of residual
    % x     = Turning points (taking values 1,...,n). [T,1]
    % n     = Number of levels.
    
    function [RFM0,res,nres] = dtp2rfm_init(n)
    
    % Initiate variables RFM0,res,nres
    
    RFM0 = zeros(n);
    res = zeros(2*n+1,1);
    nres = 0;
    
    
    function [RFM0,res,nres] = dtp2rfm1(x,RFM0,res,nres)
    
    % Count one TP.
    
    
    % Calculate RFM and res
    
    %for i = 1:length(x)
    nres = nres+1;
    res(nres) = x; %(i);
    cycleFound = 1;
    while cycleFound==1 && nres>=4
        A = sort([res(nres-1) res(nres-2)]);
        B = sort([res(nres) res(nres-3)]);
        if A(1) >= B(1) && A(2) <= B(2)
            RFM0(res(nres-2),res(nres-1)) = RFM0(res(nres-2),res(nres-1)) + 1;
            res(nres-2) = res(nres);
            nres = nres-2;
        else
            cycleFound = 0;
        end
    end
    %end
    
    
    
    function [RFM] = dtp2rfm_collect(RFM0,res,nres)
    
    % Collect RFM0 and residuals. Store in RFM.
    
    
    % Calculate RFM = RFM0 + 'cycles in res'
    
    RFM = RFM0;
    for i=1:2:nres-1
        RFM(res(i),res(i+1)) = RFM(res(i),res(i+1)) + 1;
    end
    
    % Convert to symetric rainflow
    
    RFM = RFM+RFM';
    RFM = triu(RFM);
    
