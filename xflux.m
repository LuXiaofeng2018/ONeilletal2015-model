function fa=xflux(f,u,dx,dt)
  global l r;
  fl=f(:,l);fr=f;
% include just this for centered
% also include for flux-corrected
  fa=0.5*u.*(fl+fr);
% include just this for upwind
%  fa=(v>=0).*v.*fl+(v<0).*v.*fr;
% include these for simple flux correction
%  fb=(v>=0).*v.*fl+(v<0).*v.*fr;
%  ind=(fa>0).*(fa>dx/4/dt*fl)+(fa<0).*(fa<-dx/4/dt*fr);
%  fa=(1-ind).*fa+ind.*fb;
 return