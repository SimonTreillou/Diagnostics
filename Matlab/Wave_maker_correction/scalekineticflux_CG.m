function [Ef] = scalekineticflux_CG(x,y,u,v,scales,gs)
    
    Ef=zeros(length(scales),length(y)-3,length(x)-3)+NaN;
    n=1;
    for scale_tmp=scales
        % create convolution kernel 
        radius = int16(scale_tmp/2);
        x_disk=-radius:radius+1;
        y_disk=-radius:radius+1; y_disk=y_disk';
        disk= x_disk.^2 + y_disk.^2 <= radius^2;
        disk=double(disk);
        disk=disk/(sum(sum(disk)));
    
        % convolute the fields with the kernel
        um=conv2(u,disk,'same');
        uum=conv2(u.^2,disk,'same');
        vm=conv2(v,disk,'same');
        vvm=conv2(v.^2,disk,'same');
        uvm=conv2(u.*v,disk,'same');
    
        % compute the horizontal derivatives of um and vm
        dumdx = (um(2:end-1,3:end)-um(2:end-1,1:end-2))/(gs*2);
        dumdy = (um(3:end,2:end-1)-um(1:end-2,2:end-1))/(gs*2);
        dvmdx = (vm(2:end-1,3:end)-vm(2:end-1,1:end-2))/(gs*2);
        dvmdy = (vm(3:end,2:end-1)-vm(3:end,2:end-1))/(gs*2);
    
        % compute the scale kinetic energy flux
        rho_0 = 1024.4; % define standard density in kg/m^3
        pi_tmp = -1 * rho_0 ... 
                *((uum(2:end-1,2:end-1) - um(2:end-1,2:end-1).^2).*dumdx ...
                 +(uvm(2:end-1,2:end-1) - um(2:end-1,2:end-1).*vm(2:end-1,2:end-1)).*(dumdy+dvmdx) ...
                 +(vvm(2:end-1,2:end-1) - vm(2:end-1,2:end-1).^2).*dvmdy);
        
        % set the pixels at the boundary to NaN, where the convolution kernel
        % extends over the boundary
        pi_tmp2 = zeros(size(pi_tmp)) + NaN;
        pi_tmp2(radius:end-radius,radius:end-radius) = pi_tmp(radius:end-radius,radius:end-radius);
        Ef(n,:,:)=pi_tmp2;
        n=n+1;
    end
    Ef(isnan(Ef))=0;
end