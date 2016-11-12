%KDEDEMO2 Demonstrate the difference between transformation- and ordinary-KDE
%
% KDEDEMO2 shows that the transformation KDE is a better estimate for
% Rayleigh distributed data around 0 than the ordinary KDE.
%
% Example
% kdedemo2;
%
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -*- Mode: Matlab -*- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kdedemo2.m --- 
%% Author          : Per Andreas Brodtkorb
%% Created On      : Fri Nov 19 13:29:53 2004
%% Last Modified By: Per Andreas Brodtkorb
%% Last Modified On: Sat Feb 05 10:54:25 2005
%% Update Count    : 20
%% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = rndray(1,300,1);

x = linspace(sqrt(eps),5,55);

f = kde(data);
pdfplot(f)
title('Ordinary KDE')
hold on
plot(x,pdfray(x,1),':')
hold off
 
%plotnorm((data).^(L2)) % gives a straight line => L2 = 0.5 reasonable

kopt = kdeoptset('L2',0.5);

f1 = kde(data,kopt,x);
figure(gcf+1)
pdfplot(f1)
title('Transformation KDE')
hold on
plot(x,pdfray(x,1),':')
hold off
figure(gcf-1)