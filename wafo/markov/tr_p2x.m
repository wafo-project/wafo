function [X,r] = tr_p2x(P,trans)
% TR_P2X  Transform P-matrix to X-vector
%
% Transforms a transition matrix P to a vector X containing
%   all elements in P except the diagonal of P.
%
% CALL: [X,r] = tr_p2x(P,trans)
%
% X     = Vector of length n=r*(r-1).    [nx1]
% r     = size of P-matrix.
%
% P     = transition matrix.             [rxr]
% trans = 0: No transformation. (default)
%         1: log-odds-transformation.
%              y = log(x/(1-x))
%
% See also trX2P.


if nargin<2, trans=[]; end
if isempty(trans), trans=0; end

r = length(P);
E= eye(r);
EE = E(:);
IE = find(EE==0);
PP = P';
X = PP(:);
X = X(IE);

switch trans

case 0  % No transformation

case 1  % log-odds-transformation

  X = logOdds(X);

otherwise

  error(['Transformation ' num2str(trans) ' not defined.']);

end % switch

%
% log-odds
%

function y = logOdds(x)

y=log(x./(1-x));

