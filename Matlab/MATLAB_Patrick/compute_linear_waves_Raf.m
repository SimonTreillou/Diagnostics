clear all
close all
%
%======================================================================
%  COMPUTE LINEAR SUMMATION OF WAVE SPECTRUM
%======================================================================
%
% Vyzikas et al., 2018: The evolution of free and bound waves during 
% dispersive focusing in a numerical and physical flume. Coastal 
% Engineering, 132, 95–109.
%
scale=0.154/0.05;  % scaling for flume experiment
h=1;               % Tank depth
g=9.8;
%======================================================================
%
%  Read spectrum data and compute wavenumber
%
load freq_amp_SWASH.txt

amp=freq_amp_SWASH(:,2);        % amplitude (m)
frq=freq_amp_SWASH(:,1)*2*pi;   % angular frequency (rad/s)
pha=freq_amp_SWASH(:,4);        % phase (rad)

khd=h*frq.^2/g;
kh = sqrt( khd.*khd + khd./(1.0 + khd.*(0.6666666666 ...
          +khd.*(0.3555555555 + khd.*(0.1608465608 ...
          +khd.*(0.0632098765 + khd.*(0.0217540484 ...
          +khd.*0.0065407983)))))) );
kw=kh./h;                       % wave number
%
%  Compute wave height and surface currents
%
t0=64;    % focal time
x0=14.1;  % and position

time=[0:0.1:128];

x=0;     % waves at x=0 (wave generation)
for i=1:length(time);
  ssh1(i)=sum(amp.*cos(kw*(x-x0) -frq*(time(i)-t0) -pha));
  u1(i)  =sum(amp.*frq.*cos(kw*(x-x0) -frq*(time(i)-t0) -pha))/h;
end

x=14.1;  % waves at focal point
for i=1:length(time);
  ssh2(i)=sum(amp.*cos(kw*(x-x0) -frq*(time(i)-t0) -pha));
  u2(i)  =sum(amp.*frq.*cos(kw*(x-x0) -frq*(time(i)-t0) -pha))/h;
end

%
%  Plot
%
ssh1=scale*ssh1;   % wave height
ssh2=scale*ssh2;
figure
plot(time,ssh1,time,ssh2)
legend('x=0','x=14.1')
xlabel('Time (s)')
ylabel('ssh (m)')
grid on
axis([40 70 -Inf Inf])

u1=scale*u1;       % surface currents
u2=scale*u2;
figure
plot(time,u1,time,u2)
legend('x=0','x=14.1')
xlabel('Time (s)')
ylabel('u velocity (m/s)')
grid on
axis([40 70 -Inf Inf])

figure             % spectrum
plot(frq/(2*pi),scale*amp)
xlabel('Frequency (Hz)')
ylabel('Amplitude (m)')
grid on
axis([0 1.4 -Inf Inf])






