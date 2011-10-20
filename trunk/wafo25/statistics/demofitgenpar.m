%% DEMOFITGENPAR Script to check the variance of estimated parameters
% 
% 
% Example 
% method = 'mps'
% opts.outputDir = fullfile(pwd,'html',method);
% %opts.outputDir = tempdir;
% opts.format = 'html';
% file = publish('demofitgenpar',opts);
% web(file)
%
% See also fitgenpar

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


close all
%clear all
keep method

%rand('state',sum(100*clock))  % random state each time
rand('state',1)              % same state each time

%% Parameters: Np = 30, cases = 100, k = -0.5:0.25:0.5

Np  = 30;
cases = 100;
k0 = 0.5;
k = -k0:k0/2:k0;
s = 1;
m = 0;
if ~exist('method','var'),
    method = 'mps';
end

phat = zeros(cases,3);
pstd = phat;
ciL = phat;
ciU = phat;
Nk = length(k);
phatAll = cell(1,Nk);
varAll = phatAll;

txts = {'Shape','Scale'};

for iy = 1:Nk
  h = waitbar(0,'Please wait...');
  for ix = 1:cases
    R = rndgenpar(k(iy),s,m,Np,1);
    phat0 = fitgenpar(R,'method',method);
    phat(ix,:) = phat0.params;
    pcov = phat0.covariance;
    ciL(ix,:) = phat0.lowerbound;
    ciU(ix,:) = phat0.upperbound;
    pstd(ix,:) = sqrt(diag(pcov).');
    
    waitbar(ix/cases,h)
  end
  close(h)
  phatAll{iy} = phat;
  %stdAll{iy} = pvar;
  
  
  truth1 = [k(iy),s];
  for iz = 1:2
    txt2 = txts{iz};
    headtxt = sprintf('Estimated %s (Truth = %g %g) (%s)',txts{iz},truth1,method);
    if ~(iy==1 && iz ==1)
      figure(gcf+1);clf;
    end
    r12 = phat(:,iz);
    errorbar(1:cases,phat(:,iz),2*pstd(:,iz),'b.'); hold on
    plot(1:cases,[ciL(:,iz),ciU(:,iz)],'r-')
    
    it = isfinite(phat(:,iz));
    if iz==1
      m1 = mean(phat(it,iz));
      s1 = std(phat(it,iz));
    else
      m1 = exp(mean(log(phat(it,iz))));
      s1 = std(phat(it,iz));
    end
    hline(truth1(iz),'g')
    hline([m1-2*s1, m1 m1+2*s1], 'b--'),hold off
    if 1 %k(iy)>0 || iz>1
      legend(txt2,'CI','CI','Truth',['mean of ' txt2 ' +- 2 std'],'','','Location','BestOutside')
%      legend(txt2,'CI','CI','Truth',['mean of ' txt2 ' +- 2 std'],'','','Location','SouthEast')
    else
      legend(txt2,'CI','CI','Truth',['mean of ' txt2 ' +- 2 std'])
    end
    axis([-inf,inf [-1 1]+(iz==2)+ (iz==1)*k(iy)])
    xlabel('Simulation nr')
    title(headtxt)
  end  
end

%maximizefigs

