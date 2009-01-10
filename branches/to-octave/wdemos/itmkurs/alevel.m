function L=alevel(v,l)
% ALEVEL Slice virkler.mat at a given crack length.
[n,m]=size(v);
L=[];
for j=2:m
  for i=1:n-1
    if v(i,1)<=l & v(i+1,1)>l
      L=[L (v(i,j)+v(i+1,j))/2];
    end
  end
end
