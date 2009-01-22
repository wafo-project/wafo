function R=rfadd(R1,R2)
% RFADD Addition of two rational functions
%
%   CALL:   R = rfadd(R1,R2);
%
% Addition of two rational functions. A rational function 
%
% (a(1)+a(2)*x+...a(d1)*x^(d1-1)) / (b(1)+b(2)*x+...b(d2)*x^(d2-1))
%
% should be stored as a 2 x d matrix, where d=max(d1,d2):
% [a zeros(1,d-d1);b zeros(1,d-d2)]

RR=[conv(R1(1,:),R2(2,:))+conv(R2(1,:),R1(2,:));conv(R1(2,:),R2(2,:))];
[waste,newlength]=min(find(sum(abs(RR),1)~=0));
RR(:,1:newlength-1)=[];
R=RR;