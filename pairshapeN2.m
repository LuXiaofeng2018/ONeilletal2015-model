function wlayer=pairshapeN2(locs,x,y,Br2,Wsh,N,dx)
  %solve for one gaussian shape first, with size det. by 'rad'
  %where rad is #of dx radius, probs 4
%  rad = 4;
  rad = ceil(sqrt(Br2^(-1))./dx);
  [xg,yg] = meshgrid(-rad:rad,-rad:rad);
  %xg = xg.*dt;
  %yg = yg.*dt;
  %for gaussians
  %FWHM = approx 2.355c, therefore add factor of 0.3606 to convert to Rst
  gaus=Wsh.*exp(-(Br2.*dx.^2)./0.3606*((xg+0.5).^2+(yg+0.5).^2));
  
  %create dh layer
  wlayer = zeros(size(x));
  
  buf = rad;
  bufmat = zeros(N+2.*rad,N+2.*rad);
  nlocs = locs+rad;
  
  
  %now read in center locations and apply to each center
  jj = 1;
  %but first make the center locations into corner locations
  corners = nlocs-rad;
  while jj <= size(locs,1)
  bufmat(corners(jj,1):corners(jj,1)+length(gaus)-1,...
            corners(jj,2):corners(jj,2)+length(gaus)-1)...
            = bufmat(corners(jj,1):corners(jj,1)+length(gaus)-1,...
            corners(jj,2):corners(jj,2)+length(gaus)-1)+gaus;
        jj = jj+1;
  end
  
  wlayer = bufmat(buf+1:N+buf,buf+1:N+buf);
  
  addlayer1 = zeros(size(wlayer));
  addlayer2=addlayer1;addlayer3=addlayer1;addlayer4=addlayer1;
  addcorn1=addlayer1;addcorn2=addcorn1;addcorn3=addcorn1;addcorn4=addcorn1;
  
  addlayer1(1:buf,:) = bufmat(buf+N+1:end,buf+1:N+buf);
  addlayer2(:,1:buf) = bufmat(buf+1:N+buf,buf+N+1:end);
  addlayer3(N-buf+1:N,:) = bufmat(1:buf,buf+1:N+buf);
  addlayer4(:,N-buf+1:N) = bufmat(buf+1:N+buf,1:buf);
  
  addcorn1(1:buf,1:buf) = bufmat(buf+N+1:end,buf+N+1:end);
  addcorn2(N-buf+1:end,N-buf+1:end) = bufmat(1:buf,1:buf);
  addcorn3(1:buf,N-buf+1:end) = bufmat(buf+N+1:end,1:buf);
  addcorn4(N-buf+1:end,1:buf) = bufmat(1:buf,buf+N+1:end);
  
  wlayer = wlayer+addlayer1+addlayer2+addlayer3+addlayer4...
                 +addcorn1 +addcorn2 +addcorn3 +addcorn4;
  
  %figure; pcolor(wlayer); shading flat
  layersum = sum(sum(wlayer));
  return