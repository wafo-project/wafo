function E=holmquist2(S,Q1,Q2)
% HOLMQUIST2 Computes moments for higher order reliability methods.

S=.5*(S+S');
Q1=.5*(Q1+Q1');
Q2=.5*(Q2+Q2');
E=2*trace(Q1*S*Q2*S)+trace(Q1*S)*trace(Q2*S);
