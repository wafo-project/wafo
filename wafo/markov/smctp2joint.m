function [Q,QQ] = smctp2joint(P,F)
%SMCTP2JOINT  Calculates the joint MCTP for a SMCTP.
%
% CALL: [Q,QQ] = smctp2joint(P,F)
%
% Q      = Cell array of min-max and max-min transition 
%          matrices for joint MCTP.                   {1x2}
% QQ     = Cell array of min-max and max-min transition 
%          matrices matrices for SMCTP.               {rx2}
%
% P      = Transition matrix for regime process.      [rxr]
% F      = Cell array of min-Max and Max-min matrices {rx2}
% F{i,1} = min-Max matrix, process i                  [nxn]
% F{i,2} = Max-min matrix, process i                  [nxn]
%
% If a matrix F{i,2}=[], then the process will
% be assumed to be time-reversible.
%
% See also  smctp2stat, mctp2stat

% Tested  on Matlab  5.3
%
% History:
% Updated by PJ 18-May-2000
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1999

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(2,2,ni));

% Define 

r = length(P);   % Number of regime states
n = length(F{1,1});  % Number of levels

% Check that the rowsums of P are equal to 1

P = mat2tmat(P);

% Normalize the rowsums of F{1,1},...,F{r,1} to 1
%  ==>  QQ{1,1},...,QQ{r,1}

for i = 1:r
  QQ{i,1} = F{i,1};
  QQ{i,1} = mat2tmat(QQ{i,1},1);
end

% Normalize the rowsums of F{1,2},...,F{r,2} to 1
%  ==>  QQ{1,2},...,QQ{r,2}

for i = 1:r
  
  if isempty(F{i,2})        % Time-reversible
    QQ{i,2} = F{i,1}';
  else                   % F{i,2} is given
    QQ{i,2} = F{i,2}; 
  end
    
  QQ{i,2} = mat2tmat(QQ{i,2},-1);

end

Q = cell(1,2);

% Make the transition matrix Q for the joint min-Max process

Q{1,1} = zeros(n*r,n*r);
I = 0:r:(n-1)*r;
for z=1:r
  Q0 = kron(QQ{z,1},P);
  Q{1,1}(I+z,:) = Q0(I+z,:);
end


% Make the transition matrix Qh for the joint Max-min process

Q{1,2} = zeros(n*r,n*r);
I = 0:r:(n-1)*r;
for z=1:r
  Q0 = kron(QQ{z,2},P);
  Q{1,2}(I+z,:) = Q0(I+z,:);
end
