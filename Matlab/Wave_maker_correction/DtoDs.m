function [Ds] = DtoDs(D,h)
%From depth-normalized dye concentration to surface values
%  See Hally-Rosendahl & Feddersen 2016 (Eq. 7)

hdye = min(h,2.7);
Ds = (h./hdye) .* D;
end