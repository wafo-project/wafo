function [Aa,Average]=spec2lasym(S,opt,alpha,Nsim)
%SPEC2LASYM Simulates asymmetry measures for Lagrange waves from spectrum
%
%   Useful to compute asymmetry measures for many degrees of asymmetry
%   as specified by alpha = input vector of lalpha-values
%   See help text for spec2slopedat

optl=opt;
Aa=struct('A',[],'lalpha',[],'lambda',[]);
Average=[];
Na=length(alpha);
for k=1:Na,
    optl.lalpha=alpha(k);
    Slop=spec2slopedat(S,Nsim,'time',0,optl);
    Aa.A(k)=Slop.A;
    Aa.lalpha(k)=alpha(k);
    Aa.lambda(k)=Slop.lambda;
    Average.meanwaveup{k}=Slop.meanwaveup;
    Average.meanwavedown{k}=Slop.meanwavedown;
    Average.x{k}=Slop.meanwavex;
end
