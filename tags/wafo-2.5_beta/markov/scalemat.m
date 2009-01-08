function [F1,param1] = scalemat(param,F,m,s,param1)
% SCALEMAT Scale and translate a cycle matrix.
%
% [F1,param1] = scalemat(param,F,m,s,param1)

F = flipud(F)'; % Convert to PJ-def

u=levels(param);
v=levels(param1);

n1=param1(3);
F1 = zeros(n1,n1);

for i=1:n1-1
  for j=i+1:n1
    ui = (v(i)-m)/s;
    uj = (v(j)-m)/s;
    F1(i,j)=interp2(u,u,F,uj,ui);
  end
end

[I,J] = find(isnan(F1)==1);
for k=1:length(I)
  F1(I(k),J(k)) = 0;
end

F1 = flipud(F1');
