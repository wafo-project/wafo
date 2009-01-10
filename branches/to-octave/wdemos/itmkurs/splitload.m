function [xx,x,z] = splitload(x0,tz)
% SPLITLOAD Split a wsitching load (e.g. switchingload.mat)  
% [xx,x,z] = splitload(x0,tz)
  
j0 = min(find(x0(:,1)>=tz(1,1)));
j = max(find(x0(:,1)<tz(end,1)));
x = x0(j0:j,:);
z = 3*ones(length(x),1);

j0 = 1;
for i = 1:length(tz)-1
  j1 = max(find(x(:,1)<tz(i+1,1)));
  z(j0:j1) = tz(i,2)*ones(j1-j0+1,1) ; 
  j0 = j1+1;
end
for k = min(z):max(z)
  xx{k} = x(z==k,:);
end

