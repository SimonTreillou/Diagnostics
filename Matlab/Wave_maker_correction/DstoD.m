function [D] = DstoD(Ds,h)
%From surface dye concentration to depth-normalized
%  See Hally-Rosendahl & Feddersen 2016 (Eq. 7)

hdye = min(h,2.7);
D = (hdye./h) .* Ds;
end