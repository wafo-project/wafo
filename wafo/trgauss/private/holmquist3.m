function E=holmquist3(S,Q1,Q2,Q3)
% HOLMQUIST3 Computes moments for higher order reliability methods.

S=.5*(S+S');
Q1=.5*(Q1+Q1');
Q2=.5*(Q2+Q2');
Q3=.5*(Q3+Q3');
%E1=0;
E1=8*trace(Q1*S*Q2*S*Q3*S);
E1=E1+2*trace(Q1*S*Q2*S)*trace(Q3*S)+2*trace(Q1*S*Q3*S)*trace(Q2*S)+ ...
  2*trace(Q2*S*Q3*S)*trace(Q1*S);
E=E1+trace(Q1*S)*trace(Q2*S)*trace(Q3*S);
