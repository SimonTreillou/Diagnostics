function [x] = findx_IB09(fname)
    xl=ncread(fname,'xl'); xl=xl(1);
    h=ncread(fname,'h'); h=h(:,1);
    [~,ix]=min(abs(h));
    xr=ncread(fname,'x_rho');
    xr=squeeze(xr(:,1));
    x=xr-xl+(xl-xr(ix));
end