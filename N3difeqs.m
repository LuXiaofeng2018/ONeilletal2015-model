
% %{
t=0;
tc=0;

%test
uhatvec = 0;
del2psivec = 0;
psi2vec = 0;
CFL1vec = 0;
CFL2vec = 0;


  %TIME STEPPING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if AB == 2
    u1_p=u1;v1_p=v1;h1_p=h1;
    u2_p=u2;v2_p=v2;h2_p=h2;
  end
  %---------------------------------------
  if AB == 3
    u1_p=0;v1_p=0;h1_p=0;
    u2_p=0;v2_p=0;h2_p=0;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ts=[];
hm=[]; %height max min!
psi2 = 0*x;
dhdt = psi2;
pv1 = psi2;
pv2 = psi2;
zeta1 = psi2;
zeta2 = psi2;
B2 = psi2;
B1p = B2;
pv1 = B2; pv2 = B2;

%  u  h,B
%  z  v
%

%%
%-------------------------------------------------
figure;
ii=0;
zeta1mat = [];
zeta2mat = [];
hmat = [];
Wpulsemat = [];

%}
%tpl = 25*dtinv;

while t<=tmax+dt/2;%tmax+dt/2
  
  %TIME STEPPING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if AB == 2;
    tmp=u1;u1=1.5*u1-0.5*u1_p;u1_p=tmp; %
    tmp=u2;u2=1.5*u2-0.5*u2_p;u2_p=tmp; %
    tmp=v1;v1=1.5*v1-0.5*v1_p;v1_p=tmp; %tmp holds new value, new value is updated,
    tmp=v2;v2=1.5*v2-0.5*v2_p;v2_p=tmp; %now tmp's outdated version is called old
    tmp=h1;h1=1.5*h1-0.5*h1_p;h1_p=tmp; 
    if layers == 2.5
    tmp=h2;h2=1.5*h2-0.5*h2_p;h2_p=tmp;
    end
  end
  %---------------------------------------
  if AB == 3
    if tc > 1
    u1s=u1;u1=23/12*u1-16/12*u1_p+5/12*u1_pp;u1_pp=u1_p;u1_p=u1s;
    v1s=v1;v1=23/12*v1-16/12*v1_p+5/12*v1_pp;v1_pp=v1_p;v1_p=v1s;
    h1s=h1;h1=23/12*h1-16/12*h1_p+5/12*h1_pp;h1_pp=h1_p;h1_p=h1s;
    u2s=u2;u2=23/12*u2-16/12*u2_p+5/12*u2_pp;u2_pp=u2_p;u2_p=u2s;
    v2s=v2;v2=23/12*v2-16/12*v2_p+5/12*v2_pp;v2_pp=v2_p;v2_p=v2s;
    if layers == 2.5
        h2s=h2;h2=23/12*h2-16/12*h2_p+5/12*h2_pp;h2_pp=h2_p;h2_p=h2s;
    end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % add friction
  du1dt=viscN2(u1,Re,n);
  du2dt=viscN2(u2,Re,n);
  dv1dt=viscN2(v1,Re,n);
  dv2dt=viscN2(v2,Re,n);
  
  if spongedrag1 > 0
      du1dt=du1dt-spdrag1.*(u1);
      du2dt=du2dt-spdrag2.*(u2);
      dv1dt=dv1dt-spdrag1.*(v1);
      dv2dt=dv2dt-spdrag2.*(v2);
  end
  
     % absolute vorticity
     zeta1=1-Bt.*rdist.^2+(1/dx)*(v1-v1(:,l)+u1(l,:)-u1);%used to be just 'f'
     zeta2=1-Bt.*rdist.^2+(1/dx)*(v2-v2(:,l)+u2(l,:)-u2);%
     
  % add vorticity flux, zeta*u
     zv1=zeta1.*(v1+v1(:,l));
     zv2=zeta2.*(v2+v2(:,l));

  du1dt=du1dt+0.25*(zv1+zv1(r,:));
  du2dt=du2dt+0.25*(zv2+zv2(r,:));
  
     zu1=zeta1.*(u1+u1(l,:));
     zu2=zeta2.*(u2+u2(l,:));
  
  dv1dt=dv1dt-0.25*(zu1+zu1(:,r));
  dv2dt=dv2dt-0.25*(zu2+zu2(:,r));
  
  [B1p,B2p] = BernN2(u1,v1,u2,v2,gm,c22h,c12h,h1,h2,ord);
  
  % add -grad B primes
  du1dtsq=du1dt-(1/dx)*(B1p-B1p(:,l)); %A1
  du2dtsq=du2dt-(1/dx)*(B2p-B2p(:,l)); %B1
  
  dv1dtsq=dv1dt-(1/dx)*(B1p-B1p(l,:)); %A2
  dv2dtsq=dv2dt-(1/dx)*(B2p-B2p(l,:)); %B2
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if AB == 2
  %step forward squiggle velocities
  u1sq=u1_p+dt*du1dtsq;
  u2sq=u2_p+dt*du2dtsq;
  
  v1sq=v1_p+dt*dv1dtsq;
  v2sq=v2_p+dt*dv2dtsq;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if mode == 1
      
  if rem(t,tpulseper) == 0 & t ~= 0;
      %remember this will only happen once, not for a duration
      %need to set a 'clock'
      tclock = t
      locs=paircountN2(num,N);
      wlayer=pairshapeN2(locs,x,y,Br2,Wsh,N,dx);
      newWmat=pairfieldN2(L,dx,h1,wlayer);
  end
  
  if tclock + tpulsedur > t & tclock ~= 0;
  Wmat = newWmat;
  elseif t > tpulsedur
      Wmat = 0.*x.*y;
      tclock = 0;
  end

  end
  
  Fx1=xflux(h1,u1,dx,dt)-kappa./dx.*(h1-h1(:,l));
  Fy1=yflux(h1,v1,dx,dt)-kappa./dx.*(h1-h1(l,:)); 
  dh1dt=-(1/dx)*(Fx1(:,r)-Fx1+Fy1(r,:)-Fy1);
  
  if layers == 2.5
  Fx2=xflux(h2,u2,dx,dt)-kappa./dx.*(h2-h2(:,l));
  Fy2=yflux(h2,v2,dx,dt)-kappa./dx.*(h2-h2(l,:)); 
  dh2dt=-(1./dx).*(Fx2(:,r)-Fx2+Fy2(r,:)-Fy2);
  end
  
  if tradf > 0
      dh1dt = dh1dt - 1./tradf.*(h1-1); 
      dh2dt = dh2dt - 1./tradf.*(h2-1);
  end  
%{  
  if spongedrag1 > 0
      dh1dt = dh1dt - (h1-1).*spdrag1;
  end
  if spongedrag2 > 0
      dh2dt = dh2dt - (h2-1).*spdrag2;
  end
%}
  %ADDRESS W/LAYERS, NOT SURE HOW
  if mode == 1
    dh1dt = dh1dt + Wmat;
    if layers == 2.5; dh2dt = dh2dt - H1H2.*Wmat; end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %TIME STEPPING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if AB == 3 %needs RK steps to start
    if tc <= 1 %RK steps here
    du1dt1 = u1sq+dt.*du1dtsq; %these steps already have the new
    du2dt1 = u2sq+dt.*du2dtsq; %velocity with the pressure correction
    dv1dt1 = v1sq+dt.*dv1dtsq;
    dv2dt1 = v2sq+dt.*dv2dtsq;
    dh1dt1  = h1 + dt.*dh1dt;
    if layers == 2.5; dh2dt1 = h2 + dt.*dh2dt; end
    
    u1_pp=u1_p;u1_p=u1;v1_pp=v1_p;v1_p=v1;
    u2_pp=u2_p;u2_p=u2;v2_pp=v2_p;v2_p=v2;
    h1_pp=h1_p;h1_p=h1;
    if layers == 2.5; h2_pp=h2_p;h2_p=h2; end
    
    u1sq=u1sq+dt./2.*(du1dtsq+du1dt1); %this needs to be done PRE press.corr.
    u2sq=u2sq+dt./2.*(du2dtsq+du2dt1);
    v1sq=v1sq+dt./2.*(dv1dtsq+dv1dt1);
    v2sq=v2sq+dt./2.*(dv2dtsq+dv2dt1);
    h1=h1+dt./2.*(dh1dt+dh1dt1);
    if layers == 2.5; h2=h2+dt./2.*(dh2dt+dh2dt1); end
    
    
    else
    h1=h1_p+dt.*dh1dt;
    if layers == 2.5; h2=h2_p+dt.*dh2dt; end
    end
  end
    
  if AB == 2
      h1=h1_p+dt.*dh1dt;
      if layers == 2.5; h2=h2_p+dt.*dh2dt;
         %eta = (g31.*h1+g32.*h2)./g; 
      end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      u1 = u1sq;
      u2 = u2sq;
      v1 = v1sq;
      v2 = v2sq;
  
%mach = max(max(sqrt(u1.^2+v1.^2)))./sqrt(gp*min(min(H-h)))

  %CFL1 = 4.*max(max(sqrt(u1.^2)))*dt./dx + 4.*max(max(sqrt(v1.^2)))*dt./dx;
  %CFL2 = 4.*max(max(sqrt(u2.^2)))*dt./dx + 4.*max(max(sqrt(v2.^2)))*dt./dx;
  %CFL1vec = [CFL1vec,CFL1];
  %CFL2vec = [CFL2vec,CFL2];

if rem(tc,tpl)==0
    print = 'mean h1 is '
    mean(mean(h1))
    ii=ii+1;
    ts=[ts,t];
    t
   
    u1mat(:,:,ii) = u1; v1mat(:,:,ii) = v1;
    u2mat(:,:,ii) = u2; v2mat(:,:,ii) = v2;
    h1mat(:,:,ii)  = h1;
    h2mat(:,:,ii) = h2;
    
   %filename
    %%%pcolor([zeta2./h2;zeta1./h1]); shading flat; colorbar
  %  h1mass = mean(mean(h1))
  %plot(eta(:,N/2)); hold on
  %plot(eta(:,N/2)-h1(:,N/2)); hold on
  %plot(eta(:,N/2)-h1(:,N/2)-h2(:,N/2))
  %B2phsum = mean(mean((g31-g21).*h1-g32.*h2))
  %ylim([-1.1 0.1]);
  %%%title(num2str(t));
  %%%drawnow();
  %{
  if rem(t,1000) == 0
      filehalf = strcat('T512_t_',num2str(round(t)),...
          '.mat')
     save(filehalf)
      clearvars u1mat u2mat v1mat v2mat h1mat h2mat
      ii = 0;
  end
  %}
    %if mode ==1
    %    Wpulsemat(:,:,ii)= Wmat;
%        mex1mat(:,:,ii) = mexu1 + mexv1; m1m = sum(sum(mex1mat));
%        mex2mat(:,:,ii) = mexu2 + mexv2; m2m = sum(sum(mex2mat));
    %end
    %avh = 1-mesampfreq.*dtinv;
    %mean(mean(h))
end


  %TIME STEPPING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tc=tc+1;
  t=tc*dt;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %quit program if height field breaks
  if sum(isnan(h1(:)))>0
      tossvar = 'h1 NaN error'
      return
  end
  if layers == 2.5 & sum(isnan(h2(:)))>0
      tossvar = 'h2 NaN error'
      return
  end

end
