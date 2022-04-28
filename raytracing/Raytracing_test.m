clear all
close all

dx = 10;
dz = 10;

load velocities;
vp = vel_smoother(V, 64, .05, .05, 1);
[nz, nx] = size(vp);
[gslx, gslz] = gradient(1./vp,dx,dz);

save param.mat vp gslx gslz dx dz nx nz;
xx = (0:nx-1)*dx;
zz = (0:nz-1)*dz;

sx = 5000;
sz = 3350;
vps = vp(floor(sz/dz) + 1, floor(sx/dx) + 1);
ds = 12.5;
min_angle = -80+180;
max_angle = 80+180;
na = 101;
da = (max_angle-min_angle)/na;
%%% Plot velocity model
figure, imagesc(xx,zz,vp);
T = zeros(31,1);
H = zeros(31,1);
for i=1:na
    an = min_angle + da*(i-1);
    x = [sx,sz, sind(an)/vps, cosd(an)/vps];
    %%% Integrate raytracing equations with Runge-kutta method
    [X, S] = rkn(x,ds);
    
    %%% Plot ray
    hold on, plot(X(1:4:end), X(2:4:end), 'r');
    
    %%% Calculate travel-time vs offset
    ix = floor(X(1:4:end)/dx) +1;
    iz = floor(X(2:4:end)/dz) +1;
    t = 0;
    for j=1:length(ix)
        t = t + ds/vp(iz(j),ix(j));
    end
    if(iz(end) <= 2)        
        T(i) = t;
        H(i) = X(end-3)-sx;
    end
end

%[H_sort,I]=sort(H);
figure, plot(H(T~=0),T(T~=0),'-');
ylabel('Time (s)');
xlabel('Apperture (m)');
set(gca, 'Ydir', 'reverse');



