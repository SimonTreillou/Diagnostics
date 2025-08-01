clear all
close all
 
h=7;
g=9.8;
wp=13;           % Peak period
wf=2*pi/wp;      % Frequency
wd=10.*pi/180;   % incidence angle 

khd = h*wf^2/g;  % 
kh  = sqrt(   khd*khd + khd/(1.0 + khd*(0.6666666666 ...
              +khd*(0.3555555555 + khd*(0.1608465608 ...
              +khd*(0.0632098765 + khd*(0.0217540484 ...
                                 + khd*0.0065407983)))))) );
wki=kh/h;     % wavenumber
wl=2*pi/wki;  % wavelength

L=wl/abs(sin(wd)); % domain longshore size for periodicity 

disp(['L = ',num2str(round(L))])

% To get compatible size with different angles:
%
%          an=asin(sin(a)/n)
%
% n is the factor giving smaller angles 'an' that are
% all consistent with L

a1=wd*180/pi;
a2=asin(sin(wd)/2)*180/pi;
a3=asin(sin(wd)/3)*180/pi;
disp(['compatible Y-axis Length: ', ...
num2str(round(L)),'  ',num2str(round(L*2)),'  ',num2str(round(L*3))])
disp(['compatible incidence angles: ', ...
num2str(a1),'  ',num2str(a2),'  ',num2str(a3)])


return

