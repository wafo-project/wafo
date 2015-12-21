function C=Nlevel(v,N)
% NLEVEL Slice virkler.mat at a given number of cycles.

[n,m]=size(v);

C=[];

for j=2:m
  i = min(find(v(:,j)>=N));
  if ~isempty(i)
    c = v(i,1)-(v(i,j)-N)/(v(i,j)-v(i-1,j))*(v(i,1)-v(i-1,1));
    C = [C c];
  end
end

return

