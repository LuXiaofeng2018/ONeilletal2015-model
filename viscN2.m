%visc.m
function field=viscND(vel,Re,n);

global l r l2 r2 dx

%n is exponent of Laplacian operator
%where visc term is nu * (-1)^(n+1) (\/^2)^n
% ( the \/ is my del operator)
%so for regular viscosity, n = 1
%hyperviscosity can be any n>1
%
%in this function I only permit n = 1,2

if n == 1;
    
    field = vel(l,:)+vel(r,:)+vel(:,l)+vel(:,r)-4*vel;
    field = (nu/dx^2).*field;

end

if n == 2;
    
    field = 2*vel(l,l) + 2*vel(l,r) + 2*vel(r,l) + 2*vel(r,r)...
          - 8*vel(l,:) - 8*vel(r,:) - 8*vel(:,l) - 8*vel(:,r)...
          +  vel(l2,:) +  vel(r2,:) +  vel(:,l2) +  vel(:,r2)...
          +20*vel;
    field = -1./Re.*(1/dx^4).*field;
      
end

return