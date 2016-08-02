%MOparams.m

global l r l2 r2 dx;

 %div by higher # for more plots in given time



%dtinv = 2^8;
%dt=1/dtinv;
dx=1./5;%L/N;
dt = 1./(2.^8);%dx./(10.*c1)
dtinv = 1./dt;
sampfreq = 5;
tpl = sampfreq.*dtinv;

N = ceil(L./dx);%ceil(10./9*pi*1./sqrt(Bt));
L = N.*dx;

l=[N,1:N-1]; l2=[N-1:N,1:N-2];
r=[2:N,1];   r2=[3:N,1:2];

% h B grid
[x,y]=meshgrid([0.5:N]*dx-L/2);
H = 1+0*x;
eta = 0.*x;
h1 = 0.*x+1;
h2 = 0.*x+1;

%u grid
[x,y]=meshgrid([0:N-1]*dx-L/2,[0.5:N]*dx-L/2);
u1 = 0*x.*y; u2 = u1;
%v grid
[x,y]=meshgrid([0.5:N]*dx-L/2,[0:N-1]*dx-L/2);
v1 = 0.*x.*y; v2 = v1;
% zeta grid
[x,y]=meshgrid([0:N-1]*dx-L/2,[0:N-1]*dx-L/2);
rdist = sqrt(x.^2+y.^2);
outerlim = L/2-0.5;
rlim = double(rdist<=outerlim);

sponge1 = ones(N).*max(rdist-outerlim,0);
sponge1 = sponge1./(max(max(sponge1)));
spdrag1 = spongedrag1.*sponge1;
sponge2 = ones(N).*max(rdist-outerlim,0);
sponge2 = sponge2./(max(max(sponge2)));
spdrag2 = spongedrag2.*sponge2;




