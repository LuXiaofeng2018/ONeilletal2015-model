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

%%% ND = nondimensional

c22h = 9;                           % ND 2nd baroclinic gravity wave speed squared
c12h = 10;                          % ND 1st baroclinic gravity wave speed squared
H1H2 = 1;                           % ND upper to lower layer height ratio
Bt = 1.^2./2./30.^2;                % ND scaled beta Ld2^2/4a^2
Br2 = 1.5;                          % ND scaled storm size: Burger number Ld2^2/Rst^2
p1p2 = 0.95;                        % ND upper to lower layer density ratio
tstf = 48;                          % ND storm duration tst*f0
tstpf=60;                           % ND period between forced storms tstp*f0
tradf=2000;                         % ND Newtonian damping of layer thickness trad*f0
Ar = 0.15;                          % ND areal storm coverage
Re = 5e4;                           % ND Reynolds number
Wsh = 0.03/2;                       % ND convective Rossby number

%%%%% derived quantities %%%%%

gm = p1p2.*c22h./c12h.*H1H2         % ND reduced gravity
aOLd = sqrt(1./Bt./2);              % ND planetary radius to deformation radius ratio
L = 3.*pi./9.*aOLd;         %%%???  % ND num = ceil(numfrc.*L.^2./Br2);
num = round(Ar.*L.^2.*Br2./pi)      % number of storms
deglim = 90-3.*L./2.*aOLd.*180./pi; % domain size [degrees]

iind = 1;

while iind < 2

    % use next line if running simulations sequentially from vectors above
    clearvars -except c*h H1H2 Bt Br2 p1p2 gm tstf tstpf tradf num Re Wsh L p1p2 iind

global dx;
AB = 2;
layers = 2.5;                       % or 2 for solid bottom boundary;

filename = strcat('N3_id428but_tst',num2str(tstf),...
           '_tstp',num2str(tstpf),'_dx',num2str(dx),'.mat');

spongedrag1 = 0.1;
spongedrag2 = 0.1;

n = 2;                              % order of Laplacian; '2' is hyperviscosity    
    
kappa = 1e-6;
ord = 2; %must equal 1 for Glenn's order, otherwise for Sadourney's (squares before avgs)
tmax = 8000;

N3params
N3modes
N3difeqs

%save(filename)                     % activate to save simulation output
%printthis = strcat(filename,' saved')

iind = iind+1;

end
