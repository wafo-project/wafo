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
%  F = magic(3);
%  assert(F, [ 8   1   6;...
%              3   5   7;...
%              4   9   2]);
%  assert(mat2tmat(F), ...
%         [0.5333333333333333   0.0666666666666667   0.4000000000000000;...
%          0.2000000000000000   0.3333333333333333   0.4666666666666667;...
%          0.2666666666666667   0.6000000000000000   0.1333333333333333], 1e-10);
%  assert(mat2tmat(F,1),[0.0   0.142857142857143   0.857142857142857;...
%                        0.0,  0.0, 1.0;...
%                        0.0,  0.0, 0.0], 1e-10);
%  assert(mat2tmat(F,-1), [0.0,   0.0,   0.0;...
%                          1.0,   0.0,   0.0;...
%                          0.307692307692308   0.692307692307692   0.0], 1e-10);
%  assert(mat2tmat(F,-1,-2), [ 0   0   0; 0   0   0; 1   0   0], 1e-10);
%  assert(mat2tmat(F,-1,0),... 
%                  [1.0,   0.0,   0.0;...
%                   0.3750,   0.6250,   0.0;...
%                   0.266666666666667,   0.60,   0.133333333333333], 1e-10);

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (P�r Johannesson) 1998
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.1'

% Check input arguments

ni = nargin;
no = nargout;
%error(nargchk(1,3,ni));
narginchk(1,3)
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
