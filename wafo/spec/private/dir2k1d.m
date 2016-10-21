function Snew=dir2k1d(S)
  
  Snew=S;
  g=gravity;
  % Commence with adding S(.,theta)+S(.,theta+pi) for -pi/2<theta<pi/2
  % since instantaneously there is no difference between waves coming
  % from direction -pi/4 and 3*pi/4 for example.
  dth=diff(S.theta(1:2));
  th0=min(S.theta(S.theta>=0));
  th=[fliplr(th0-dth:-dth:-pi/2) th0:dth:pi/2]';
  [r2,c]=size(S.S(S.theta>pi/2+1e-13,:));
  [r3,c]=size(S.S(S.theta<-pi/2-1e-13,:));

  m1=length(th);
  if r2==0 % no thetas over pi/2, which is the largest theta?
    m1=find(abs(S.theta(end)-th)<=1e-8);
  end
  m2=1;
  if r3==0 % no thetas under -pi/2, which is the smallest theta?
    m2=find(abs(S.theta(1)-th)<=1e-8);
  end
  S0=zeros(length(th),length(S.w));
  S0(m2:m1,:)=S.S(r3+1:end-r2,:);
  S0(1:r2,:)=S0(1:r2,:)+S.S(end-r2+1:end,:);
  S0(end-r3+1:end,:)=S0(end-r3+1:end,:)+S.S(1:r3,:);

  % Divide the integral in two parts, abs(th)<=a0 and abs(th)>a0 
  k=linspace(0,S.w(end)^2/g,length(S.w))';
  a0=pi/4;
  % the case abs(th)<=a0
  th1=th(th>=-a0 &th<=a0);
  S1=S0(th>=-a0 &th<=a0,:);
  k1=k(2:end,ones(1,length(th1)));
  th1=th1(:,ones(1,length(S.w)-1))';
  wmax=max(max(sqrt(g*k1./cos(th1))));
  S.w=[S.w; S.w(end)+diff(S.w(end-1:end));wmax*1.1];
  S1=[S1 zeros(size(S1,1),2)]; % add zeros to prevent NaN in interpolation
  Sip=interp2(S.w,th1(1,:),S1,sqrt(g*k1./cos(th1)),th1);
  Snew.S=[0; sqrt(g./k(2:end)).*simpson(th1(1,:),Sip'./sqrt(cos(th1))',1).'/2];
  
  % the case abs(th)>a0
  th1=th;
  S1=S0;
  wmax=sqrt(g*k(end)/cos(a0));
  n=2;
  if wmax*1.1>S.w(end)
    S.w=[S.w; wmax*1.1];
    n=3;
  end
  S1=[S1 zeros(size(S1,1),n)]; % add zeros to prevent NaN in interpolation
  S1=[S1(end,:); S1; S1(1,:)]; % extend angles to prevent NaN in interpolation
  th1=[th1(1)-dth;th1;th1(end)+dth];
  Sk=zeros(size(k));
  for j=2:length(k)
    w=S.w(S.w>=sqrt(g*k(j)/cos(a0)));
    thw=acos(g*k(j)./w.^2);
    Sip=(interp2(S.w,th1,S1,w,thw)+(interp2(S.w,th1,S1,w,-thw)))./sqrt(w.^4-(g*k(j))^2);
    if length(w(Sip>0))==2
      Sk(j)=mean(Sip(Sip>0))*diff(w(Sip>0));
    elseif length(w(Sip>0))>2
      Sk(j)=g*simpson(w(Sip>0),Sip(Sip>0),1);
% not perfect here, somtimes change of integration limit in simpson
    end
  end
  hold off
  plot(k,Snew.S,k,Sk)
 
  Snew.S=Snew.S+Sk;
  Snew=rmfield(Snew,'w');
  Snew=rmfield(Snew,'theta');
  Snew.k=k;
  Snew.type='rotk1d';
  spec2mom(Snew) 
  hold on
  plotspec(Snew,1,'r')