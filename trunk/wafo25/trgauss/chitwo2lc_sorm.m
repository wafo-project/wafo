function  mu = chitwo2lc_sorm(u,g,b,S22,S12,ass,n0,epsi)
%CHITWO2LC_SORM  SORM-approximation of crossing intensity, noncentral Chi^2 process 
%
%  The noncentral Chi^2 process is defined in IR-(27).
%           X(t)=sum beta_j Z_j(t) + gamma_j Z_j(t)^2 
%
%  CALL:  mu = chitwo2lc_sorm(x0,gamma,beta,S12,S22,ass);
%  
%       mu   = two column vector containing: levels in the first column, 
%              the SORM-approximation of crossing intensitity for quadratic sea  
%              in the second column.
%       u    = column with levels, note that u can not take value zero.
%     g,b    = the vectors containing constants defining the process
%              approximation of the cdf of the quadratic sea.  
%     S12    = covariance between vector Z(t) and Z'(t), where Z(t)=(Z_1(t),...,Z_n(t)).
%     S22    = covariance between vector Z'(t) and Z'(t)., where Z(t)=(Z_1(t),...,Z_n(t)).
%     ass    = 0, SORM, 0< FORM, >0 higher order very slow,  optional ass=0;
%      n0    = parameter used to find a starting point in the optimization, default value n0=10
%    epsi    = tolerance level in the optimization, default value epsi=1e-9
%
% Example: 
%   S=jonswap; 
%   [gamma beta S12 S22]=dirsp2chitwo(S.S,S.w);
%   L0=sqrt(sum(beta.^2)); u=(0.01:0.1:4*L0)'; u=[(-4*L0:0.1:-0.1)';u];
%   mu = chitwo2lc_sorm(u,gamma,beta,S22,S12); 
%   semilogy(mu(:,1),mu(:,2))
%
% See also chitwo2lc_sp, dirsp2chitwo


% References: 
%             K. Breitung, (1988), "Asymptotic crossing rates for stationary {G}aussian vector
%                                   processes.", Stochastic Processes and
%                                   their Applications, 29,pp. 195-207.
%                                   
%             U. Machado, I. Rychlik (2002) "Wave statistics in nonlinear sea", Extremes, 6, pp. 125--146.
%             Hagberg, O. (2005) PhD - thesis, Dept. of Math. Statistics, Univ. of Lund. 
% Baxevani, A.,  Hagberg, O. and Rychlik, I.  Note on the distribution of extreme waves crests (OMAE 2005).


%   Calls: mindist
% Revised by I.R 14.04.05
% By OH.         24.10.04
%
%------------------------------------------------------------------------------------

if nargin<7
   n0=10;
end 
if nargin<8
   epsi=1e-9;
end
%
% if ass>0 one is using second order assymptotics, see PhD thesis of Hagberg
%
if nargin<6
   ass=0;
end

n=length(u);
x0 = zeros(length(g),n);
%x0=[];
mu=[];
for i=1:n;
    xx1=mindist(g,b,u(i),n0,epsi);
    x0(:,i)= xx1;
    if ass>0
        [Cb0 Cb1]=rqlf_asympt(u(i),g,b,S22,S12);
        %mu=[mu; u(i) muu];
        mu=[mu; u(i) 0.5*(Cb0+Cb1)];
    end
end
if ass<=0
    
g=g(:);
b=b(:);
d=length(g);
[d1,n]=size(x0);
if d1~=d
   error('the observations of x0 should be stored as columns')
end
b0=sqrt(sum(x0.^2,1));

G=g*(-2*b0./sqrt(sum((kron(b,ones(1,n))+2*kron(g,ones(1,n)).* ...
		      x0).^2)));
% cofactor matrix:
D=zeros(size(x0));
for k=1:d
   D(k,:)=prod(1+G([1:k-1 k+1:d],:),1);
end

%cc=[sqrt((sum(x0.*(S22*x0),1)+sum(G.*(S12'*x0).^2,1))) sqrt(sum(D.*x0.^2,1))];

mu=exp(-b0.^2/2).*sqrt((sum(x0.*(S22*x0),1)+sum(G.*(S12'*x0).^2,1))./(sum(D.*x0.^2,1)))/pi/2;

%mu1=exp(-b0.^2/2);%.*sqrt((sum(x0.*((S22)*x0),1))./(sum(x0.^2,1)))/2/pi;
mu2=exp(-b0.^2/2).*sqrt((sum(x0.*((S22)*x0),1))./(sum(x0.^2,1)))/2/pi;

%cc=sqrt((sum(x0.*(S22*x0),1)+sum(G.*(S12'*x0).^2,1))./(sum(D.*x0.^2,1)))/pi/2;
mu=[u mu'];
if ass<0
 mu=[u mu2'];
end
end