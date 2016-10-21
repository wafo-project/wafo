function y=cdfmargcnd2dfun(H,V)
%CDFMARGCND2DFUN is an internal function to cdfmargcnd2d prbmargcnd2d. 
%
%  CALL: y = cdfmargcnd2dfun(x2,x1)
%
% Dependending on condon it returns a product of conditional 
% pdf and cdf.
%   condon = 0  returns p(X1,X2)=p(X2)*P( X1|X2) 
%            1  returns p(X1)*p( X1|X2) 
%            2  returns P( X1|X2) 
%            3  returns p( X1|X2) 
%
%
%   X1 and X2 must have equal size.
%   The size of P is the common size of the arguments X1 and X2.
%   
% GLOBALS used PHAT  CONDON
%
% See also   cdfmargcnd2d, prbmargcnd2d

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


% tested on: matlab 5.2
% history:
% pab 09.11.99


global PHAT CONDON
UDIST=lower(PHAT.dist{2});
CDIST=lower(PHAT.dist{1});
PH=PHAT.x{2};


[Av , Bv,Cv]=margcnd2dsmfun(PHAT,H);
switch CONDON
  case {0,1} , % no conditional or conditional CDF given V
   switch UDIST(1:2)
    case 'ra',  pdf1= pdfray(H,PH);
    case 'we' ,  pdf1=pdfweib(H,PH(1),PH(2));
    case 'gu' ,  pdf1=pdfgumb(H,PH(1),PH(2),0);
    case 'tg' ,  pdf1=pdfgumb(H,PH(1),PH(2),1);
    case 'ga' ,  pdf1=pdfgam(H,PH(1),PH(2));
    case 'gg',  pdf1=pdfgengam(H,PH(1),PH(2),PH(3));
    case 'lo' ,  pdf1=pdflognorm(H,PH(1),PH(2));
    otherwise, error('unknown distribution')
   end 
  case {2,3}, pdf1=1;%conditional CDF given H
end

switch CONDON
  case {0,2}
    switch CDIST(1:2)
      case 'ra', y=pdf1.*cdfray(V-Cv,Av);
      case 'gu',y =  pdf1.*cdfgumb(V-Cv,Av,Bv,0);
      case 'tg',  y =  pdf1.*cdfgumb(V-Cv,Av,Bv,1);
      case 'lo', y =  pdf1.*cdflognorm(V-Cv,Av,Bv);
     case 'ga', y =  pdf1.*cdfgam(V-Cv,Av,Bv);	
      case 'gg', y =  pdf1.*cdfgengam(V,Av,Bv,Cv);	
      case 'we', y =  pdf1.*cdfweib(V-Cv,Av,Bv);
      otherwise, error('Unknown distribution')
    end
  case {1,3},
   switch CDIST(1:2)
    case 'ra', y=pdf1.*pdfray(V-Cv,Av);
    case 'gu',y =  pdf1.*pdfgumb(V-Cv,Av,Bv,0);
    case 'tg',  y =  pdf1.*pdfgumb(V-Cv,Av,Bv,1);
    case 'lo', y =  pdf1.*pdflognorm(V-Cv,Av,Bv);
    case 'ga', y =  pdf1.*pdfgam(V-Cv,Av,Bv);	
    case 'gg', y =  pdf1.*pdfgengam(V,Av,Bv,Cv);
    case 'we', y =  pdf1.*pdfweib(V-Cv,Av,Bv);
    otherwise, error('Unknown distribution')
   end
end

