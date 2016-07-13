%N2energies.m

ii = 1;%size(u1mat,3);
kj = 1;
while ii < size(u1mat,3)+1

    h1n = h1mat(:,:,ii);
    h2n = h2mat(:,:,ii);
    
K1t(kj,1) = sum(sum(1./2.*gm.*c12h./c22h.*Axl(Ayl(h1n)).*(u1mat(:,:,ii).^2+v1mat(:,:,ii).^2))).*dx.^2./L.^2;
K2t(kj,1) = sum(sum(1./2.*Axl(Ayl(h2n)).*(u2mat(:,:,ii).^2+v2mat(:,:,ii).^2))).*dx.^2./L.^2;

At(kj,1) = sum(sum(1./2.*gm.*c12h./c22h.*c12h.*(Axl(Ayl(h1n))-1).^2....
+1./2.*c22h.*(Axl(Ayl(h2n))-1).^2+gm.*c12h.*(Axl(Ayl(h1n))-1).*(Axl(Ayl(h2n))-1))).*dx.^2./L.^2;

ii = ii + 1;
kj = kj + 1;

end
%figure
%plot(ts,K1t+K2t,ts,At,ts,K1t+K2t+At,'k')


p1p2 = gm.*c12h./c22h./H1H2;
Ar = num.*pi./L.^2./Br2;
Ast = (0.5*H1H2.*c22h+0.5*p1p2.*c12h+gm.*c12h).*H1H2....
       *(Wsh.*tstf).^2.*(Ar./(1-Ar));
   
Ep = Ast.*tradf./tstpf;
Etimefrac = round(0.8*length(ts));
Etot = mean(K1t(Etimefrac:end)+K2t(Etimefrac:end)+At(Etimefrac:end));
Atot = mean(At(Etimefrac:end));

%clearvars u1mat u2mat v1mat v2mat h1mat h2mat
