function plotspec(S,varargin),
%PLOTSPEC Plot a spectral density 
%
% CALL:  plotspec(S,plotflag,linetype) 
%
%    S       = an array of spectral density structs:
%              1D  (see dat2spec)
%              2D  (see createspec)
%   1D:
%   plotflag = 1 plots the density, S, (default)
%              2 plot 10log10(S)
%	       3 plots both the above plots 
%   2D:
%   Directional spectra: S(w,theta), S(f,theta)             
%   plotflag = 1 polar plot S (default)
%              2 plots spectral density and the directional 
%                spreading, int S(w,theta) dw or int S(f,theta) df
%              3 plots spectral density and the directional 
%                spreading, int S(w,theta)/S(w) dw or int S(f,theta)/S(f) df
%              4 mesh of S
%              5 mesh of S in polar coordinates  
%              6 contour plot of S
%              7 filled contour plot of S
%   Wavenumber spectra: S(k1,k2)
%   plotflag = 1 contour plot of S (default)
%              2 filled contour plot of S
%   lintype : specify color and lintype, see PLOT for possibilities.
%
% NOTE: - lintype may be given anywhere after S.
%
% Examples
%  S = demospec('dir'); 
%  S2 = mkdspec(jonswap, spreading);
%  plotspec(S,2); hold on;
%  plotspec(S,3,'g');  % Same as previous fig. due to frequency independent spreading
%  plotspec(S2,2,'r'); % Not the same as previous figs. due to frequency dependent spreading
%  plotspec(S2,3,'m');
%  % transform from angular frequency and radians to frequency and degrees
%  Sf = ttspec(S,'f','d'); clf
%  plotspec(Sf,2);
%
%  close all;
%
% See also  dat2spec, createspec, simpson

% Changed reference to  wfindpeaks  (GL Feb 2011)
% tested on Matlab 5.3
% History:
% revised pab 2008
% -renamed from wspecplot to plotspec
% Revised pab feb2007
% -Fixed bug: S.phi was interpreted as counter-clockwise rotation, but is
%             clockwise rotation
% Revised pab oct2005
% -Fixed bug: S.phi was interpreted as clockwise rotation, but is
%            counter-clockwise rotation
% Revised pab Apr2005
%  Fixed a bug when S is a multidimensional struct.
% revised pab Feb2004
%  - fixed abug in fcolorbar  
% revised pab dec2003
% replaced call to colorbarf with fcolorbar 
% replaced call to findzlevel with clevels  
% revised pab 29.06.2001
% directional spectra:
% - added filled contour plot (made also for wave number spectra)
% - clarified what is actually plotted for plotflag==2 and plotflag==3
%   + added examples on the subject
% - ylabels for plotflag==5 and S.theta given in radians fixed.
% revised pab 22.06.2001
%  -added contourplot for directional spectra + examples in help header.
%  -Fixed bugs: S.phi was not taken into acount for plotflag 2:5
%               and also phi was not transformed to radians for
%               plotflag=1 when spectrum was given in degrees.
% revised by IR 02.04.2001: l 191. Before plot not worked when Fp=0.
% revised by pab 15.01.2001
%  - fixed a bug: axis for plotflag==2
% revised by gl 12.07.2000
%  - added third argument 0 on line 125 to prevent warning from findrfc
% revised pab 14.06.2000
% - added the possibility that S.theta is given in degrees
% - fixed bug for plotflag 2 S.type = dir. (simpson integration)
% - cleaned up some code and put it into separate functions:
%        num2pistr, findzlevel cltext1 and fixthetalabels
% revised pab 08.05.2000
%  - fixed a bug for S an array of structs
% revised jr 10.04.2000
%  - made a comment of line 127
% revised pab 28.01.2000
%  - added plotflag 4 for dir. spec. 2D
%  - added point-and-click editing of all the text objects (title, xlabel, ylabel) of the current figure
% revised es 21.01.2000 lintype also for dir.spec.   
% revised pab 20-01.2000
%  - added the possibility that S is a array of structs
% revised pab 12-01.2000
%  - added plotflag 3 for dir. spec. i.e., meshplot
%  - made the xticks and xticklabels print better 
%    for dir. spec. when plotflag==2 (line 244-257)
% revised pab 06.01.2000
%  - added checking for empty c from contourc line 190
% revised pab 22.10.1999
%  - fixed a bug: legend(txt) -> legend(txt{:}) line 147  
% revised pab 12.10.1999
%
% by pab,es 17.09.1999

% label the contour levels
txtFlag = 0;
LegendOn = 1;

P   = varargin;

ih=ishold;
Ns=length(S);
if Ns>1
  cfig=gcf;
  for ix=1:Ns,
    if ih
      newplot;
    else
      figure(cfig-1+ix);
    end
    plotspec(S(ix),P{:});
  end
  return
end

% Default values: 
%~~~~~~~~~~~~~~~
plotflag = 1;
lintype  = 'b-';
%Nlevel   = 7;

Np    = length(P);
if Np>0
  strix = zeros(1,Np);
  for ix=1:Np, % finding symbol strings 
    strix(ix)=ischar(P{ix});
  end
  k = find(strix);
  if any(k)
    Nk = length(k);
    if Nk>1
      warning('WAFO:plotspec','More than 1 strings are not allowed');
    end
    lintype = P{k(1)};
    Np = Np-Nk;
    P  = {P{find(~strix)}}; % remove strings from input
  end
  if Np>0 && ~isempty(P{1}), plotflag=P{1};end
end


ftype = freqtype(S); %options are 'f' and 'w' and 'k'
switch ftype
case {'w'},
  freq      = S.w; 
  xlbl_txt  = 'Frequency [rad/s]';  
  ylbl1_txt = 'S(w) [m^2 s / rad]';
  ylbl3_txt = 'Directional Spectrum';
  zlbl_txt  = 'S(w,\theta) [m^2 s / rad^2]';
  funit     = ' [rad/s]';  
  %Sunit     = ' [m^2 s / rad]'; 
case {'f'}, 
  freq      = S.f; 
  xlbl_txt  = 'Frequency [Hz]';     
  ylbl1_txt = 'S(f) [m^2 s]';
  ylbl3_txt = 'Directional Spectrum';
  zlbl_txt  = 'S(f,\theta) [m^2 s / rad]';
  funit     = ' [Hz]';  
  %Sunit     = ' [m^2 s ]'; 
case {'k'},
  freq      = S.k; 
  xlbl_txt  = 'Wave number [rad/m]';
  ylbl1_txt = 'S(k) [m^3/ rad]';
  funit     = ' [rad/m]'; 
  %Sunit     = ' [m^3 / rad]';
  if isfield(S,'k2')
    ylbl4_txt='Wave Number Spectrum';
  end
otherwise,  error('Frequency type unknown');
end

if isfield(S,'norm') && S.norm 
  %Sunit=[];
  funit=[];
  ylbl1_txt = 'Normalized Spectral density';
  ylbl3_txt = 'Normalized Directional Spectrum';
  ylbl4_txt = 'Normalized Wave Number Spectrum';
  if strcmp(ftype,'k')
    xlbl_txt = 'Normalized Wave number';
  else
    xlbl_txt = 'Normalized Frequency'; 
  end
end

ylbl2_txt = 'Power spectrum (dB)';

phi =0;
if isfield(S,'phi') && ~isempty(S.phi)
  phi = S.phi;
end


switch lower(S.type(end-2:end))
  case {'enc','req','k1d'} % 1D plot
    S.S  = S.S(:);
    Fn   = freq(end); % Nyquist frequency
    indm = wfindpeaks(S.S,4);
    if isempty(indm)
      %disp('indm empty')  % Commented, 00.04.10 
      [maxS,indm] = max(S.S);  %to be removed when wfindpeaks works
    end
    maxS = max(S.S(indm));
    if isfield(S,'CI') && ~isempty(S.CI),
      maxS  = maxS*S.CI(2);
      txtCI = [num2str(100*S.p), '% CI'];
    end
    
    Fp = freq(indm); %peak frequency/wave number
    
    if length(indm)==1
      txt = {['fp = ' num2str(Fp,2), funit];};
    else
      txt = cell(1,length(indm));
      txt(1) = {['fp1 = ' num2str(Fp(1),2), funit]};
      for j=2:length(indm)
         txt(j) = {['fp',num2str(j),' = ', num2str(Fp(j),2), funit]};
      end
    end
    if (plotflag == 3), subplot(2,1,1); end 
    if (plotflag == 1) || (plotflag ==3),% Plot in normal scale
      plot([Fp(:) Fp(:)]',[zeros(length(indm),1) S.S(indm)]',':');
      hold on;
      if isfield(S,'CI'),
        plot(freq,S.S*S.CI(1), 'r:' );
        plot(freq,S.S*S.CI(2), 'r:' );
      end
      plot(freq,S.S,lintype);
      if ~ih, hold off;end
      if ih, a=axis; else a=zeros(1,4); end
      a1=Fn;
      if (Fp>0) 
        a1=max(min(Fn,10*max(Fp)),a(2));
      end
      axis([0 a1 0 max(1.01*maxS,a(4))]);
      title('Spectral density');
      xlabel(xlbl_txt);
      ylabel(ylbl1_txt);
    end
    
    if (plotflag==3), subplot(2,1,2); end
    
    if (plotflag == 2) || (plotflag ==3), % Plot in logaritmic scale
      ind=find(S.S>0);
      
      plot([Fp(:),Fp(:)]',[repmat(min(10*log10(S.S(ind)/maxS)),length(Fp),1) 10*log10(S.S(indm)/maxS)]',':');
      hold on;
      if isfield(S,'CI'),
        plot(freq(ind),10*log10(S.S(ind)*S.CI(1)/maxS), 'r:' );
        plot(freq(ind),10*log10(S.S(ind)*S.CI(2)/maxS), 'r:' );
      end
      plot(freq(ind),10*log10(S.S(ind)/maxS),lintype);
      if ~ih, hold off; end
      if ih, a=axis; else a=[0 0 0 0]; end
      axis([0 max(min(Fn,max(10*Fp)),a(2)) -20 max(1.01*10*log10(1),a(4))]); % log10(maxS)
      title('Spectral density');
      xlabel(xlbl_txt);
      ylabel(ylbl2_txt);
    end      
    if LegendOn
      if isfield(S,'CI'),
        legend(txt{:},txtCI);
      else
        legend(txt{:});
      end
    end
  case {'k2d'}
    if plotflag==1,
      [c, h] = contour(freq,S.k2,S.S,'b');
      z_level = clevels(c);

      if txtFlag==1
        textstart_x=0.05; textstart_y=0.94;
        cltext1(z_level,textstart_x,textstart_y);
      else
        cltext(z_level,0);
      end
    else
      [c,h] = contourf(freq,S.k2,S.S);
      %clabel(c,h), colorbar(c,h)
      fcolorbar(c); % alternative
    end
    rotate(h,[0 0 1],-phi*180/pi);

    xlabel(xlbl_txt);
    ylabel(xlbl_txt);
    title(ylbl4_txt);
    %return
    km = max([-freq(1) freq(end) S.k2(1) -S.k2(end)]);
    axis([-km km -km km]);
    hold on;
    plot([0 0],[ -km km],':');
    plot([-km km],[0 0],':');
    axis('square');
    
    
    %cltext(z_level);
    %axis('square')
    if ~ih, hold off; end
  case {'dir'}
    thmin = S.theta(1)-phi;
    thmax = S.theta(end)-phi;
    if plotflag==1, % polar plot
      if 0, % alternative but then z_level must be chosen beforehand
        h = polar([0 2*pi],[0 freq(end)]);
        delete(h);hold on;
        [X,Y] = meshgrid(S.theta,freq);
        [X,Y] = pol2cart(X,Y);
        contour(X,Y,S.S',lintype);
      else
        if (abs(thmax-thmin)<3*pi), % angle given in radians
          theta = S.theta;
        else
          theta = S.theta*pi/180; % convert to radians
          phi  = phi*pi/180;
        end
        c = contourc(theta,freq,S.S');%,Nlevel); % calculate levels
        if isempty(c)
          c = contourc(theta,freq,S.S);%,Nlevel); % calculate levels
        end
        [z_level c] = clevels(c); % find contour levels
        h = polar(c(1,:),c(2,:),lintype);
        rotate(h,[0 0 1],-phi*180/pi);
      end
      title(ylbl3_txt);
      % label the contour levels
      
      if txtFlag==1,
        textstart_x = -0.1; textstart_y=1.00;
        cltext1(z_level,textstart_x,textstart_y);
      else
        cltext(z_level,0);
      end
      
    elseif (plotflag==2) || (plotflag==3),
      %ih = ishold;
      
      subplot(211);
      
      if ih, hold on; end
      
      Sf = spec2spec(S,'freq'); % frequency spectrum
      plotspec(Sf,1,lintype);

      subplot(212);
      
      Dtf        = S.S;
      [Nt,Nf]    = size(S.S); 
      Sf         = Sf.S(:).';
      ind        = find(Sf);
      
      if plotflag==3, %Directional distribution  D(theta,freq))
        Dtf(:,ind) = Dtf(:,ind)./Sf(ones(Nt,1),ind);
      end
      Dtheta  = simpson(freq,Dtf,2); %Directional spreading, D(theta)
      Dtheta  = Dtheta/simpson(S.theta,Dtheta); % make sure int D(theta)dtheta = 1
      [y,ind] = max(Dtheta);
      Wdir    = S.theta(ind)-phi; % main wave direction
      txtwdir = ['\theta_p=' num2pistr(Wdir,3)]; % convert to text string
      
      plot([1 1]*S.theta(ind)-phi,[0 Dtheta(ind)],':'); hold on;
      if LegendOn,
        lh=legend(txtwdir,0);
      end
      plot(S.theta-phi,Dtheta,lintype);
      
      fixthetalabels(thmin,thmax,'x',2);  % fix xticklabel and xlabel for theta
      ylabel('D(\theta)');
      title('Spreading function');
      if ~ih, hold off; end
      %legend(lh) % refresh current legend
    elseif plotflag==4, % mesh
      mesh(freq,S.theta-phi,S.S);
      xlabel(xlbl_txt);
      fixthetalabels(thmin,thmax,'y',3); % fix yticklabel and ylabel for theta
      zlabel(zlbl_txt);
      title(ylbl3_txt);
    elseif plotflag==5, % mesh
      %h=polar([0 2*pi],[0 freq(end)]);
      %delete(h);hold on
      [X,Y]=meshgrid(S.theta-phi,freq);
      [X,Y]=pol2cart(X,Y);
      mesh(X,Y,S.S');
      % display the unit circle beneath the surface
      hold on; mesh(X,Y,zeros(size(S.S')));hold off;
      zlabel(zlbl_txt);
      title(ylbl3_txt);
      set(gca,'xticklabel','','yticklabel','');
      lighting phong;
      %lighting gouraud
      %light
    elseif (plotflag==6) || (plotflag==7),
      theta = S.theta-phi;
      [c, h] = contour(freq,theta,S.S); %,Nlevel); % calculate levels
      fixthetalabels(thmin,thmax,'y',2); % fix yticklabel and ylabel for theta
      if plotflag==7,
        hold on;
        [c,h] =	contourf(freq,theta,S.S); %,Nlevel); % calculate levels
        %hold on
      end
    
      title(ylbl3_txt);
      xlabel(xlbl_txt);
      if 0,
        [z_level] = clevels(c); % find contour levels
        % label the contour levels
        if txtFlag==1
          textstart_x = 0.06; textstart_y=0.94;
          cltext1(z_level,textstart_x,textstart_y); % a local variant of cltext
        else
          cltext(z_level);
        end
      else
        colormap('jet');
        if plotflag==7,
          fcolorbar(c);
        else
          %clabel(c,h),
          hcb = colorbar;
        end
        grid on;
      end
    else
      error('Unknown plot option');
    end
  otherwise, error('unknown spectral type');
end

if ~ih, hold off; end

%  The following two commands install point-and-click editing of
%   all the text objects (title, xlabel, ylabel) of the current figure:

%set(findall(gcf,'type','text'),'buttondownfcn','edtext')
%set(gcf,'windowbuttondownfcn','edtext(''hide'')')

return

function  xtxt = num2pistr(x,N)
% NUM2PISTR convert a scalar x to a text string in fractions of pi
%           if the numerator is less than 10 and not equal 0 
%           and if the denominator is less than 10.
%
%  CALL xtxt = num2pistr(x,n)
% 
% xtxt = a text string in fractions of pi
%  x   = a scalar
%  n   = maximum digits of precision. (default 3)


if nargin<2||isempty(N)
  N=3;
end
den = 0;
num = 0;
if x~=0,
  [num, den] = rat(x/pi);
end
if (den<10) && (num<10) && (num~=0),
  if abs(den)==1,
    dtxt=''; 
  else
    dtxt=['/' num2str(den)];
  end % denominator
  if abs(num)==1, % numerator
    if num==-1,
      ntxt='-';
    else
      ntxt='';
    end
  else 
    ntxt=num2str(num); 
  end
  xtxt= [ntxt '\pi' dtxt];
else
  xtxt = num2str(x,N);
end
return  
  
function fixthetalabels(thmin,thmax,xy,dim)
%FIXTHETALABELS pretty prints the ticklabels and x or y labels for theta  
%  
% CALL fixthetalabels(thmin,thmax,xy,dim)
%
%  thmin, thmax = minimum and maximum value for theta (wave direction)
%  xy           = 'x' if theta is plotted on the x-axis 
%                 'y' if theta is plotted on the y-axis 
%  dim          = specifies the dimension of the plot (ie number of axes shown 2 or 3)
%  If abs(thmax-thmin)<3*pi it is assumed that theta is given in radians 
%  otherwise degrees

ind = [('x' == xy),  ('y' == xy) ];
yx = 'yx';
yx = yx(ind);
if nargin<4||isempty(dim),
  dim=2;
end
%drawnow
%pause

if abs(thmax-thmin)<3*pi, %Radians given. Want xticks given as fractions  of pi
  %Trick to update the axis 
  if xy=='x'	
    if dim<3,
      axis([thmin,thmax 0 inf ]);
    else
      axis([thmin,thmax 0 inf 0 inf]);
    end
  else
    if dim<3, 
      axis([0 inf thmin,thmax ]);
    else
      axis([0 inf thmin,thmax 0 inf]);
    end
  end
  
  set(gca,[xy 'tick'],pi*(thmin/pi:0.25:thmax/pi));
  set(gca,[xy 'ticklabel'],[]);
  x = get(gca,[xy 'tick']);
  y = get(gca,[yx 'tick']);
  y1 = y(1);
  dy = y(2)-y1;
  yN = y(end)+dy;
  ylim = [y1 yN];
  dy1 = diff(ylim)/40;
  %ylim=get(gca,[yx 'lim'])%,ylim=ylim(2);
  
  if xy=='x',
    for j=1:length(x),
      xtxt = num2pistr(x(j));
      figtext(x(j),y1-dy1,xtxt,'data','data','center','top');
    end
   % ax = [thmin thmax 0 inf];
    ax = [thmin thmax ylim];
    if dim<3,
      figtext(mean(x),y1-7*dy1,'Wave directions (rad)','data','data','center','top');
    else
      ax = [ax  0 inf];      
      xlabel('Wave directions (rad)');
    end
  else
    %ax = [0 inf thmin thmax];
    ax = [ylim thmin thmax];
   
    if dim<3,
      for j=1:length(x)
        xtxt = num2pistr(x(j));
        figtext(y1-dy1/2,x(j),xtxt,'data','data','right');
      end
      set(gca,'DefaultTextRotation',90);
      %ylabel('Wave directions (rad)')
      figtext(y1-3*dy1,mean(x),'Wave directions (rad)','data','data','center','bottom');
      set(gca,'DefaultTextRotation',0);
    else
      for j=1:length(x),
        xtxt = num2pistr(x(j));
        figtext(y1-3*dy1,x(j),xtxt,'data','data','right');
      end
      ax = [ax 0 inf];
      ylabel('Wave directions (rad)');
    end
  end
  %xtxt = num2pistr(x(j));
  %for j=2:length(x)
  %  xtxt = strvcat(xtxt,num2pistr(x(j)));
  %end
  %set(gca,[xy 'ticklabel'],xtxt)
else % Degrees given
  set(gca,[xy 'tick'],thmin:45:thmax);
  if xy=='x',
    ax=[thmin thmax 0 inf];
    if dim>=3,      ax=[ax 0 inf];    end
    xlabel('Wave directions (deg)');
  else
    ax=[0 inf thmin thmax ];
    if dim>=3,      ax=[ax 0 inf];    end
    ylabel('Wave directions (deg)');
  end
end
axis(ax);
return

function cltext1(z_level,textstart_x,textstart_y)
%CLTEXT1 Places contour level text in the current window
%
%  [h,ax] = cltext1(levels,x,y)
%
%         x       = the starting x-coordinate of the text.
%         y       = the starting y-coordinate of the text.
%

%  textstart_x = -0.07; textstart_y=1.00;

  delta_y     = 1/33;
  %dir:
  %h = figtext(textstart_x-0.07,textstart_y,' Level curves at:','normalized');

  % k2d:
  h = figtext(textstart_x-0.05,textstart_y,' Level curves at:','normalized');
  set(h,'FontWeight','Bold');

  textstart_y = textstart_y-delta_y;
  for ix=1:length(z_level)
    textstart_y = textstart_y-delta_y;
    figtext(textstart_x,textstart_y,num2str(z_level(ix),4),'normalized');
  end
return


