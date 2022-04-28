function [ X,S ] = rkn(x, h)

load param.mat

%RKN Rundge-Kutta method for ray tracing
s = 0;
S = s;
X = x';
while(1)
    k1 = gradfun(s,x, vp, nx, nz, dx, dz, gslx, gslz);
    if(sum(isnan(k1)))
        break;
    end
    k2 = gradfun(s+h/2,x+h*k1/2, vp, nx, nz, dx, dz, gslx, gslz);
    if(sum(isnan(k2)))
        break;
    end
    k3 = gradfun(s+h/2,x+h*k2/2, vp, nx, nz, dx, dz, gslx, gslz);
    if(sum(isnan(k3)))
        break;
    end
    k4 = gradfun(s+h,x+h*k3, vp, nx, nz, dx, dz, gslx, gslz);
    if(sum(isnan(k4)))
        break;
    end
    k = (k1+2*k2+2*k3+k4)/6;
    
    s = s + h;
    x = x + h*k;
    
    S = [S; s];
    X = [X; x'];
    
end

end

