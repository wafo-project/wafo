function P = tr_x2p(X,trans)
% TR_X2P  Transform X-vector to transition matrix P.
%
% CALL: P = tr_x2p(X,trans)
%
% P     = transition matrix.           [rxr]
% 
% X     = Vector of length n=r*(r-1).  [nx1]
% trans = 0: No transformation. (default)
%         1: inverse log-odds-transformation.
%              x = 1/(exp(-y)+1)
%
% See also tr_p2x

if nargin<2, trans=[]; end
if isempty(trans), trans=0; end

r=(1+sqrt(1+4*length(X)))/2;

switch trans

case 0  % No transformation

case 1  % inverse-log-odds-transformation

  X = invlogOdds(X);

otherwise

  error(['Transformation ' num2str(trans) ' not defined.']);

end % switch

E= eye(r);
EE = E(:);
IE = find(EE==0);
Y = zeros(r*r,1);
Y(IE) = X;
P = reshape(Y,r,r)';
for i = 1:r
  P(i,i) = 1-sum(P(i,:));
end

%
% inverse of log-odds
%

function x=invlogOdds(y)

x=1./(exp(-y)+1);

