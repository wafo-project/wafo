function P = mat2tmat(F,def,K)
% MAT2TMAT  Converts a matrix to a transition matrix.
%
% CALL: P = mat2tmat(F)
%       P = mat2tmat(F,def)
%       P = mat2tmat(F,def,Nsub)
%
% P   = Transition matrix			[nxn]
%
% F   = Matrix					[nxn]
% def =  0: Markov chain transition matrix. (default)
%        1: min-Max transition matrix.
%           (Zeros on and below the diagonal.)
%       -1: Max-min transition matrix.
%           (Zeros on and above the diagonal.)
% K   = 'if def=1' : Set zeros below the K-th diagonal.
%       'if def=-1': Set zeros above the K-th diagonal.
%       'if def=0' : Not applicable.
%
% The routine converts a matrix to a transition matrix, 
% i.e. normalizes each row sum to 1.
%
% Example: 
%   F = magic(5)
%   mat2tmat(F)
%   mat2tmat(F,1)
%   mat2tmat(F,-1)
%   mat2tmat(F,-1,-3)
%   mat2tmat(F,-1,0)

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1998
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.1'

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,3,ni));

if ni == 1
  def = 0;
end
if ni < 3
  K = def;
end

n = length(F);  % Number of levels

if def == 0  % Transition matrix for Markov chain
  
  P = F;
  % Set negative elements to zero
  [I,J] = find(P<0);
  if length(I) ~= 0
    warning(['Negative elements in F. Setting to zero!']);
    for k = 1:length(I)
      P(I(k),J(k)) = 0;
    end
  end

  sumP = sum(P,2);
  % Normalize rowsums
  if sum(sumP == 1) ~= n
    for i = 1:n
      if sumP(i)~=0
	P(i,:) = P(i,:)/sumP(i); 
      end
    end
  end
  
elseif def == 1 % min-Max transition matrix 
  
  P = triu(F,K); % Set zeros on and below diagonal
  % Set negative elements to zero
  [I,J] = find(P<0);
  if length(I) ~= 0
    warning(['Negative elements in F. Setting to zero!']);
    for k = 1:length(I)
      P(I(k),J(k)) = 0;
    end
  end

  sumP = sum(P,2);
  % Normalize rowsums
  N = min([n-K n]); % Number of sums that should be equal to 1.
  if sum(sumP == 1) ~= N
    for i = 1:N
      if sumP(i)~=0
	P(i,:) = P(i,:)/sumP(i); 
      end
    end
  end
  
elseif def == -1 % Max-min transition matrix 
  
  P = tril(F,K); % Set zeros on and above diagonal
  % Set negative elements to zero
  [I,J] = find(P<0);
  if length(I) ~= 0
    warning(['Negative elements in F. Setting to zero!']);
    for k = 1:length(I)
      P(I(k),J(k)) = 0;
    end
  end

  sumP = sum(P,2);
  N = min([n+K n]);
  if sum(sumP == 1) ~= N
    for i = n-N+1:n
      if sumP(i)~=0
        P(i,:) = P(i,:)/sumP(i); 
      end
    end
  end
  
else
  
  error(['Undefined definition: def = ' num2str(def) ]);
  
end
