

clearvars u*mat v*mat h*mat zeta*mat *tsq ii ts tmp tclock t tc
clearvars sampfreq tmax

sampfreq = 20;
tmax = 8000;

%This current version of the model has:
% 2 and 2.5 model option
% hyperviscosity
% radiative relaxation
% 2nd order AB time stepping
% sponge drag
% baroclinic drag

global dx;
AB = 2;
layers = 2.5; %or 2;



filename2 = strcat(filename,'_rerun.mat');

n = 2; %order of Laplacian; '2' is hyperviscosity    

ord = 2; %must equal 1 for Glenn's order, otherwise for Sadourney's (squares before avgs)


N3paramsrerun
N3modesrerun
N3difeqsrerun

save(filename2)
printthis = strcat(filename,' saved')




