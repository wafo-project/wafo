function E=holmberg1(S,a,b,Q1)
% HOLMBERG1 Computes moments for higher order reliability methods.
%
% CALL: E=holmberg1(S,a,b,Q1);
%
% Computation of the expectation of 
% sqrt(pi/2)*(b'*X)*(X'*Q*X)*(2*cdfnorm(a'*X)-1)
% if S is normally distributed with mean zero and covariance matrix
% S.
Q1=(Q1+Q1')/2;
term1=a'*S*b*trace(S*Q1);
term2=2*a'*S*Q1*S*b;
term3=-(a'*S*b)*(a'*S*Q1*S*a)/(1+a'*S*a);
E=(term1+term2+term3)/sqrt(1+a'*S*a);
