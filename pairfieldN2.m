function Wmat=pairfieldN2(L,dx,h1,wlayer)
  
voldw = sum(sum(wlayer)).*dx.^2;
area = L.^2;
wcorrect = voldw./area;
Wmat = wlayer-wcorrect;
%figure; pcolor(hlayer); hold on; shading flat

  return
