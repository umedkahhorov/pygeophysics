function [ xp ] = gradfun(s,x, vp, nx, nz, dx, dz, gslx, gslz)
%GRADFUN Function to compute dx/ds and dp/ds
%   

xp = x;
ix = floor(x(1)/dx) + 1;
iz = floor(x(2)/dz) + 1;
if(iz < 1  || iz > nz)
    xp(:) = nan;
    return;
end
if(ix < 1 || ix > nx )
    xp(:) = nan;
    return;
end
xp(1) = vp(iz,ix).*x(3);
xp(2) = vp(iz,ix).*x(4);
xp(3) = gslx(iz,ix);
xp(4) = gslz(iz,ix);
end

