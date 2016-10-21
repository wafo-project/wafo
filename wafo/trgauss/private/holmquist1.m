function E=holmquist1(S,Q1)
% HOLMQUIST1 Computes moments for higher order reliability methods.

S=.5*(S+S');
Q1=.5*(Q1+Q1');
E=trace(Q1*S);
