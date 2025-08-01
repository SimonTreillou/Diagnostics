%%Windows for Signal Processing
%%Ara Anbarasu
%function [out] = window(win_size,type)
%Input                                  |%Output
%win_size=   no:of points in window     |%out  =  window of length [win_size] 
%type    =   window type                 
%Example Usage
%win = window(2^4,'rectangular')
%Available Windows
%blackman
%flattop
%hamming
%nuttall
%rectangular

function [out] = window(win_size,type)

res = ones(win_size,1);

switch lower(type)

    case 'blackman'

        for i=0:1:win_size-1

            res(i+1) = 0.42 - 0.5*cos(2*pi*i/(win_size-1)) + 0.08*cos(4*pi*i/(win_size -1));

        end

    case 'flattop'

        a0=1;

        a1=1.93;

        a2=1.29;

        a3=0.338;

        a4=0.032;

        for i=0:1:win_size-1

            res(i+1) = a0 - a1*cos(2*pi*i/(win_size-1)) + a2*cos(4*pi*i/(win_size -1))- a3*cos(6*pi*i/(win_size -1)) + a4*cos(8*pi*i/(win_size -1));

        end

    case 'hanning'

        for i=0:1:win_size-1

            res(i+1) = 0.5*(1 - cos(2*pi*i/(win_size-1)));

        end

    case 'hamming' 

          for i=0:1:win_size-1

            res(i+1) = 0.54 - 0.46*cos(2*pi*i/(win_size-1));

          end

    case 'nuttall' 

        a0=0.3635819;

        a1=0.4891775;

        a2=0.1365995;

        a3=0.0105411;

        a4=0.032;

        for i=0:1:win_size-1

            res(i+1) = a0 - a1*cos(2*pi*i/(win_size-1)) + a2*cos(4*pi*i/(win_size -1))- a3*cos(6*pi*i/(win_size -1));

        end

    otherwise

end

out = res;

end

