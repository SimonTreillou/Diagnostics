function res = PSD(x)
% MMA.m
% Author: Alexander Khaustov; alexander dot khaustov at gmail dot com 
% Copyright (C) 2008 St.-Petersburg Institute of Cardiological Technics (Incart), www.incart.ru
% This software is released under the terms of the GNU General
% Public License (http://www.gnu.org/copyleft/gpl.html).
% 
% A candidate for proper periodogram normalization

len = length(x);
f = fft(x); 
ac = real(ifft(f .* conj(f)/ len)) ;    % autocorrelation
xx = abs(fft(ac .* hamming(len)));
%xx = abs(fft(ac));
res = xx(1:(len / 2 + 1)) / len; 

return;
