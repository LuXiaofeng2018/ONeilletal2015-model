locs=paircountN2(num,N);
mode = 1;

switch(mode)
    
    case 1 %mass flux W field
      pulse = 'off';%'on'; %'off'
      %clear locs; locs = round([L./2,3.*L./4;...%single, control
      %                    3.*L./4,L./2; 3.*L./4+Ld2,L./2;...
      %                    L./2,L./4;    L./2,L./4-1.5;...
      %                    L./4,L./2;    L./4-2.*sqrt(Br2.^(-1)),L./2]./dx)
      
      
      wlayer=pairshapeN2(locs,x,y,Br2,Wsh,N,dx);
      Wmat=pairfieldN2(L,dx,h1,wlayer);%eta-h1
      Wmatorig = Wmat;
      tpulseper = tstpf;%50; %'s', not dt
      tpulsedur = tstf; %'s', not dt
      tclock = 0;
      FreeGrid = length(find(spdrag1==0))./(N.^2);
      % h B grid
      [x,y]=meshgrid([0.5:N]*dx-L/2);
      H = 1+0*x;
      if layers == 2;
           h = 0.*x+0.5;
           h10 = 0.5;
      else
           eta = 0.*x;
           %h10 = 0.2;%-FreeGrid.*mean(mean(Wmat)).*tpulsedur./tpulseper.*tradrel1-0.005;
           h1 = 0.*x+1; 
           %h20 = 0.2;%+FreeGrid.*mean(mean(Wmat)).*tpulsedur./tpulseper.*tradrel1+0.005;
           h2 = 0.*x+1;
      end 

      %u grid
      [x,y]=meshgrid([0:N-1]*dx-L/2,[0.5:N]*dx-L/2);
      %rdist = sqrt(x.^2+y.^2);
      %rlim = double(rdist<=3);
      %%%u1 = 0.*x.*y; u2 = u1; 
      %v grid
      [x,y]=meshgrid([0.5:N]*dx-L/2,[0:N-1]*dx-L/2);
      %rdist = sqrt(x.^2+y.^2);
      %rlim = double(rdist<=3);
      %%%v1 = 0.*x.*y;
      %0.1*sin(2*pi*x./L);%0.*x.*y; 
      %%%v2 = 0.*x.*y;%v1;
      
    case 2 %balanced pairs on interface
        %fmat=0.*fmat+0.1;
      locsh = [0,-8]+N/2;
      locsc = [0,8]+N/2;
      wlayerhot=pairshape(locsh,x,y,rad,amp,N);
      Wmathot=-pairfield(L,h1,dx,wlayerhot);
      wlayercold = pairshape(locsc,x,y,rad,amp,N);
      Wmatcold = pairfield(L,h1,dx,wlayercold);
      h1 = -Wmathot-Wmatcold+0.5;
      horig = h1;
      %u grid
      %[x,y]=meshgrid([0:N-1]*dx-L/2,[0.5:N]*dx-L/2);
      finv = f.^(-1);
      u1 = -gp.*(1/dx).*finv.*(h1-h1(l,:));
      %imagesc(u1)
      u2 = -u1;%-
      %v grid
      %[x,y]=meshgrid([0.5:N]*dx-L/2,[0:N-1]*dx-L/2);
      v1 = gp.*(1/dx).*finv.*(h1-h1(:,l));
      v2 = -v1;%-

    case 3 %invert a PV field
              
      %solve for LaPlacian(psi2)
      wlayer  = pairshape(locs,x,y,rad,amp);
      htop = h-0.2*pairfield(L,h,dx,wlayer);
      h = H-htop;
      %figure; pcolor(h); hold on; shading flat
      del2psit = h;
             %else simply put in my desired pv field

      % calculate psi2 by inverting LaPlacian(psi2)
      psitop = ifft2(-1./wv2nonz.*fft2(del2psit));
      %psi2 = ifft2(-wv2inv*fft2(del2psi2));
      %hbot = H-htop(shift,:)+0.5;
      %del2psib = (hbot);
      %psibot = ifft2(-1./wv2nonz.*fft2(del2psib));
      
      u1 =  1/dx*(psitop-psitop(l,:));
      %u1 = u1(shift,:);
      v1 = -1/dx*(psitop-psitop(:,l));%-
      %v1 = v1(shift,:);
      u2 = 1/dx*(psitop-psitop(l,:));%-
      v2 =  -1/dx*(psitop-psitop(:,l));
      %figure; quiver(x,y,u1,v1,4)
      
    case 6 %field of geostrophically balanced like-signed peaks
        
      wlayer=pairshape(locs,x,y,rad,amp);
      Wmat=-pairfield(L,h,dx,wlayer);
      h = -Wmat+0.5;
      horig = h;
      %u grid
      %[x,y]=meshgrid([0:N-1]*dx-L/2,[0.5:N]*dx-L/2);
      u1 = -gp./f*(1/dx)*(h-h(l,:));
      u2 = -u1;%-
      %v grid
      %[x,y]=meshgrid([0.5:N]*dx-L/2,[0:N-1]*dx-L/2);
      v1 = gp./f*(1/dx)*(h-h(:,l));
      v2 = -v1;%-
      
    case 7 %test for linearity
        
      [x,y]=meshgrid([0.5:N]*dx-L/2);
      H=1+0.*x;
      %waveL = 8
      waveL = waveLnum*dx;
      kw = 2*pi/waveL; %wavenumber
      h1 = 0.05*sin(kw.*y)+h1;
      u1 = 0.*h1; v1 = u1; u2 = u1; v2 = u1;
      
      
    case 8 %velocity field, no height perturbation
        
      H = 1+0.*x;
      [x,y]=meshgrid([0:N-1]*dx-L/2,[0.5:N]*dx-L/2);
      u1 = 0.05+0.*x;
      %u1(N/2+1:end,:) = -0.005*x(N/2+1:end,:);
      u2 = -u1;
      v1 = 0.*x; v2 = v1;
      
    otherwise
      
end


%%
% if pulse == 'on'
%     T = 15;
%     tconv = 3;
%     pulsevec = zeros(T,1);
%     pulsevec(1:tconv,1) = ones(1:tconv,1);
%     j = 1;
%     dtlength = 0;
%     pulsesignal = zeros(T*dtinv,1);
%     while j < T+1
%         longpulsevec = repmat(pulsevec(j),1,dtinv);
%         pulsesignal(dtlength+1:dtlength+dtinv) = longpulsevec;
%         dtlength = dtlength+dtinv;
%         j = j+1;
%     end
% end
        
    
        
        
        
    
    
