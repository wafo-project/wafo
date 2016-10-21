function [l,res] = disper2(t,dt,niter,dx_rel)
%DISPER2 Dispersion relation with possible mean flow
%
%     disper2: dispersion relation, Newton-Raphson
%              with possible mean flow
%
%              (2.pi.f-k.v)^2 = g.k.tanh(k.d)
%
%     syntax:  [l,res] = disper2(t[,d[,niter[,dx_rel]]])
%
%     input:   t  : period vector
%              dt : water depth vector [water depth d, mean flow v] 
%                   or matrix same number of lines as length of t
%                   vector or matrix [water depth d, surface current ...
%                             speed, bottom current speed, type of profile]
%                   type of profile = 0 => uniform (bottom current ... 
%                             speed is not used) (default value)
%                   type of profile = 1 => linear
%                   type of profile = 2 => exponential
%                   (default value [Inf,0])
%
%       output:     l : wavelength (same size as t)
%                   res : residue (om0-xk*v)^2-g0*xk*tanh(xk*d)

% WAVEMOD Project - June 1994 - (Agnes Robin)
%                  (transfered from macro basile AR -> matlab)
% modified (mean flow) Marc Prevosto - October 2007
% modified (shear flow) Marc Prevosto - June 2015
%

vt = 0 ; v_bottomt = 0 ; type_profilet = 0 ;
if nargin < 2 || isempty(dt), dt = Inf ; end
if nargin < 3, niter = 50 ; end
if nargin < 4, dx_rel = 0.001 ; end
if min(size(dt))==1, dt = reshape(dt,1,length(dt)) ; end
if size(dt,2)==2, vt = dt(:,2) ; dt = dt(:,1) ;
elseif size(dt,2)==3, vt = dt(:,2) ; v_bottomt = dt(:,3) ;
elseif size(dt,2)==4, vt = dt(:,2) ; v_bottomt = dt(:,3) ; type_profilet = dt(:,4) ; dt = dt(:,1) ;
end

m = length(t) ;
g0 = 9.81 ;
l = zeros(m,1) ; res = zeros(m,1) ;

for iv=1:m
    ti = t(iv) ;
    if length(dt)==1, d = dt(1) ; else d = dt(iv) ; end
    if length(vt)==1, v = vt(1) ; else v = vt(iv) ; end
    if ti==0
        xl = 0 ; ff = 0 ;
    elseif isinf(ti), xl = Inf ; ff = 0 ;
    elseif isinf(d)
        if v==0
            xl = g0*ti^2/2/pi ;
        else
            om0 = 2*pi/ti ;
            discr = g0*(g0+4*om0*v) ;
            k = (g0+2*om0*v-sqrt(discr))/(2*v^2) ; xl = 2*pi/k ;
        end
        ff = 0 ;
    else
        fi = 1/ti ;
        if   fi==0
            xl = NaN ;
        else
            % first iteration with v = 0
            om0 = 2*pi*fi ;
            k0 = 4.0243*fi^2 ;
            xk = k0 ;
            ftest = 99 ;
            for ii = 1:niter,
                z   = xk*d ;
                y   = tanh(z) ;
                ff  = om0^2-g0*xk*y ;
                dff = g0*(z*(y^2-1)-y) ;
                xk_old = xk ;
                xk  = xk_old-ff/dff ;
                ftest = abs((xk-xk_old)/xk_old) ;
                if ftest <= dx_rel, break ; end
            end
            % second iteration with v uniform
            if v~=0
                ftest = 99 ;
                for ii = 1:niter
                    z = xk*d ;
                    y = tanh(z) ; om0mxkv = om0-xk*v ;
                    ff = om0mxkv^2-g0*xk*y ;
                    dff = g0*(z*(y^2-1)-y)-2*v*om0mxkv ;
                    xk_old = xk ;
                    xk = xk_old-ff/dff ;
                    ftest = abs((xk-xk_old)/xk_old) ;
                    if ftest <= dx_rel, break ; end
                end
            end
            if length(type_profilet)==1, type_profile = type_profilet(1) ; else type_profile = type_profilet(iv) ; end
            if type_profile==1
                if length(v_bottomt)==1, v_bottom = v_bottomt(1) ; else v_bottom = v_bottomt(iv) ; end
                ftest = 99 ;
                for ii = 1:niter
                    z = xk*d ; y = tanh(z) ; om0mxkv = om0-xk*v ;
                    ff = om0mxkv*(om0mxkv-(v_bottom-v)/d*y)-g0*xk*y ;
                    dff = g0*(z*(y^2-1)-y)-v*2*om0mxkv+(v_bottom-v)*(v*y/d-om0mxkv*(1-y^2)) ;
                    xk_old = xk ;
                    xk = xk_old-ff/dff ;
                    ftest = abs((xk-xk_old)/xk_old) ;
                    if ftest <= dx_rel, break ; end
                end
            elseif type_profile==2
                if length(v_bottomt)==1, v_bottom = v_bottomt(1) ; else v_bottom = v_bottomt(iv) ; end
                ftest = 99 ;
                for ii = 1:niter
                    z = xk*d ; om0mxkv = om0-xk*v ;
                    Ke = -log((om0-v_bottom*xk)/om0mxkv)/z ; sKe = sqrt(Ke^2+1) ;
                    dKe = om0*(v_bottom-v)/(om0-v_bottom*xk)/om0mxkv/d/xk-Ke/xk ;
                    y = tanh(sKe*z) ; dy = (z*Ke*dKe/sKe+sKe*d)*(1-y^2) ;
                    ff = om0mxkv^2*(sKe-Ke*y)-g0*xk*y ;
                    dff = -2*om0mxkv*(sKe-Ke*y)*v+om0mxkv^2*(Ke*dKe/sKe-dKe*y-Ke*dy)-g0*(y+xk*dy) ;
                    xk_old = xk ;
                    xk = xk_old-ff/dff ;
                    ftest = abs((xk-xk_old)/xk_old) ;
                    if ftest <= dx_rel, break ; end
                end
            else
            end
            xl = 2*pi/xk ;
            %       fprintf(' ------> %d iterations, dk_rel = %6.4f, f_k = %5.2e  \n', [ii ftest abs(ff)])
        end
    end
    l(iv) = xl ; res(iv) = ff ;
end
