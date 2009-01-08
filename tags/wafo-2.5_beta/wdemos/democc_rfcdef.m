function ccrfc=democc_rfcdef(proc,tp,point,ccrfc)
%DEMOCC_RFCDEF illustrates the definition of rainflow cycles.
%
% CALL: ccrfc=democc_rfcdef(proc,tp,point,ccrfc)
%
% Used by democc 
  
% Tested on: matlab 5.3
% History:
% Revised by PJ 10-Jan-2000
%   updated for WAFO
% Original version by Mats Frendahl
  
ms = 20; % markersize
  
democc_markmax(proc,tp,point,0);

proc=[(1:length(proc))' proc(:,2)];
time=proc(:,1); L=proc(:,2); n=length(L);

level=tp(point,2); refpoint=tp(point,1);
hold on, title('Definition of rainflow count')
crossleft=max(find(L(1:tp(point,1))>tp(point,2)));
if isempty(crossleft),
   crossleft=1;
   noleft=1;
else
   noleft=0;
end
leftdepth=min(L(crossleft:refpoint));
c1=[crossleft refpoint]; c2=[level level];
plot(c1,c2,'--'), 

pause(1)

leftindex=find(L(crossleft:refpoint)==leftdepth)+crossleft-1;
plot(leftindex,leftdepth,'k.','markersize',ms)
plot([leftindex leftindex],[leftdepth level],'--');

pause(1)

list=find(L>=level);
if ( (length(list)>1) & ( refpoint<max(list)))
   index=find(list==refpoint);
   crossright=list(index+1);
   noright=0;
else
   crossright=n;
   noright=1;
end

plot([tp(point,1) crossright],[level level],'--'), 

pause(1)

if noright == 1
   rightdepth = -Inf;
else
   rightdepth=min(L(refpoint:crossright));
   rightindex=find(L(refpoint:crossright)==rightdepth)+refpoint-1;
   plot(rightindex,rightdepth,'k.','markersize',ms)
   plot([rightindex rightindex],[rightdepth level],'--');
end

pause(1)

if leftdepth>rightdepth
   plot([leftindex leftindex],[leftdepth level]);
   rfcmin=leftdepth;
   if noright~=1,
      plot(rightindex,rightdepth,'k.','erase','xor','markersize',ms)
   end
else
   plot([rightindex rightindex],[rightdepth level]);
   rfcmin=rightdepth;
%   if noleft~=1
      plot(leftindex,leftdepth,'k.','erase','xor','markersize',ms)
%   end
end

hold off

if isempty(ccrfc)
   ccrfc=[rfcmin level];
else
   ccrfc=[ccrfc; rfcmin level];
end

