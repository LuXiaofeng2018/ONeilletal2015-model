function [B1,B2] = BernN2(u1,v1,u2,v2,gm,c22h,c12h,h1,h2,ord)
  global l r;
  if ord == 1 %Glenn's way
  B1 = 'broke';% (g31.*h1+g32.*h2) + 0.125*(u1+u1(:,r)).^2+0.125*(v1+v1(r,:)).^2;
  B2 = 'broke';%(g31-g21).*h1+g32.*h2 + 0.125*(u2+u2(:,r)).^2+0.125*(v2+v2(r,:)).^2;
  else        %Sadourney's way
  B1 = c12h.*h1+c22h.*h2+0.25*(u1.^2+u1(:,r).^2+v1.^2+v1(r,:).^2);
  
  B2 = gm.*c12h.*h1+c22h.*h2+0.25*(u1.^2+u1(:,r).^2+v1.^2+v1(r,:).^2);
  
  end
  
 return