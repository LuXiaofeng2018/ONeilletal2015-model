%ONrunmodel
clear all
%!cd /Users/morganoneill/Documents/Saturn/GlennSWCode/swper_MO

%This current version of the model has:
% 2 and 2.5 model option
% hyperviscosity
% radiative relaxation
% 2nd order AB time stepping
% sponge drag
% baroclinic drag

c22h = 9;%-0.6135;
c12h = 10;%4.044;%2.0283;
H1H2 = 1;
Bt = 1.^2./2./30.^2;           %scaled beta Ld2^2/4a^2
Br2 = 1.5;%2^2./(1.^2);%Burger Ld2^2/Rst^2
p1p2 = 0.95;      %nondim density ratio
tstf = 48;                   %tst*f0
tstpf=60;%6.25;
tradf=2000;
Ar = 0.15;
Re = 5e4;
Wsh = 0.03/2;
%%%%%
gm = p1p2.*c22h./c12h.*H1H2
aOLd = sqrt(1./Bt./2);
L = 3.*pi./9.*aOLd; %num = ceil(numfrc.*L.^2./Br2);
num = round(Ar.*L.^2.*Br2./pi)
deglim = 90-3.*L./2.*aOLd.*180./pi;

iind = 1;

while iind < 2

    clearvars -except c*h H1H2 Bt Br2 p1p2 gm tstf tstpf tradf num Re Wsh L p1p2 iind

global dx;
AB = 2;
layers = 2.5; %or 2;


filename = strcat('N3_id428but_tst',num2str(tstf),...
           '_tstp',num2str(tstpf),'_dx',num2str(dx),'.mat');

spongedrag1 = 0.1;%0.1;
spongedrag2 = 0.1;%0.1;

n = 2; %order of Laplacian; '2' is hyperviscosity    
    
kappa = 1e-6;
ord = 2; %must equal 1 for Glenn's order, otherwise for Sadourney's (squares before avgs)
tmax = 8000;

N3paramsA
N3modesA
N3difeqsA

save(filename)
printthis = strcat(filename,' saved')
iind = iind+1;

end

