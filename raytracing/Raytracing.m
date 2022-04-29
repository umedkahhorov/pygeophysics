clear all
[v,h] = ReadRSS('../vp3_global.rss');
v1 = squeeze(v);
vp = v1';
%%
dx = 10;
dz = 10;

l1_a = linspace(-40,-14,4);
l2_a = linspace(-13,-8,401);
mid_a = linspace(-8.9,8.9,4);
r1_a = linspace(8,13,401);
r2_a = linspace(14,40,5);

n1=4;
n2=401;
angs = zeros(1,n1*3+n2*2);
angs(1,1:n1) = l1_a;
angs(1,n1:n2+n1-1) = l2_a;
angs(1,n2+n1:n2+n1*2-1) = mid_a;
angs(1,n2+n1*2:n2*2+n1*2-1) = r1_a;
angs(1,n2*2+n1*2:n1*3+n2*2) = r2_a;
na = n1*3+n2*2;

[nz, nx] = size(vp);
dx = 20;
dz = 20;
[gslx, gslz] = gradient(1./vp,dx,dz);

save param.mat vp gslx gslz dx dz nx nz;
xx = (0:nx-1)*dx;
zz = (0:nz-1)*dz;

sx = 290000;
sz = 10;
vps = vp(floor(sz/dz) + 1, floor(sx/dx) + 1);
ds = 12.5;


% min_angle = 20;
% max_angle = 40;
% da = (max_angle-min_angle)/na;
% na = 21;

%%% Plot velocity model
figure, imagesc(xx,zz,vp);
T = zeros(31,1);
H = zeros(31,1);
results = zeros(na,6000);
for i=1:na
%     an = min_angle + da*(i-1);
    an = angs(1,i);
    x = [sx,sz, sind(an)/vps, cosd(an)/vps];
    %%% Integrate raytracing equations with Runge-kutta method
    [X, S] = rkn(x,ds);
  
    idx = size(X,1);
    results(i,1:idx) = X;
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



