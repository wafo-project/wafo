function [ro,PP]=mc2stat(P)
% MC2STAT  Calculates the stationary distribution for a Markov chain.
%
% CALL: ro = mc2stat(P);
%
% ro = Stationary distribution.   [1xr]
%
% P  = transition matrix.         [rxr]
%
% Example: 
%   F = magic(5)
%   P = mat2tmat(F)
%   ro = mc2stat(P)
%
% See also  mat2tmat, mc2reverse, smc2stat, mctp2stat
  
% Tested  on Matlab  5.3
%
% History:
% Revised by PJ  23-Nov-1999
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1998
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.0'
%   previous name: 'statf_mc'

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(1,1,ni));

n=length(P); % number of states

PP=P;

% Set negative elements to zero
[I,J] = find(PP<0);
if length(I) ~= 0
  warning(['Warning: Negative elements in P. Setting to zero!']);
  for k = 1:length(I)
    PP(I(k),J(k)) = 0;
  end
end

J =[];
NoDeleted = 0; NoDeleted0 = -1; 
while NoDeleted ~= NoDeleted0

  % Check that the rowsums of PP are equal to 1

  sumPP = sum(PP');
  if sum(abs(sumPP-1) <= n*eps) ~= length(PP)  % Number of rowsums equal to 1
    warning(['Warning: Rowsums of P not equal to 1. Renormalizing!']);
    for i = 1:length(PP)
      if sumPP(i) == 0
        if isempty( find(J==i) )  % Not allready deleted
          warning(['Warning: Rowsum ' num2str(i) ' of P equals zero. Deleting state ' num2str(i) '!']);
          J = [J i];
          PP(:,i) = zeros(n,1);
        end
      else
        PP(i,:) = PP(i,:)/sumPP(i);
      end
    end
  end

  NoDeleted0 = NoDeleted;
  NoDeleted = length(J);

end

I=1:n; I(J)=[]; 
P1 = PP(I,I);  % Delete states


% Compute stationary distribution

r=size(P1,1);
P1=P1';
M=[eye(r-1) zeros(r-1,1)];
P1=P1(1:r-1,:);
P1=P1-M;
P1=[P1 ; ones(1,r)];
a=[zeros(r-1,1); 1];
roP1=P1\a;
%roP1=pinv(P1)*a;
roP1=roP1';

% Set ro(i)=0 for 'singular states'

ro = zeros(1,n);
ro(I)=roP1;

