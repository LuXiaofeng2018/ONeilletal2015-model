function locs=paircountN2(num,N)
%r = rand([m,n]) returns an m-by-n matrix
%let n=2 for coordinate pair
  %locs=round(rand(num,2)*(N-2*rad-1)+rad+1);
  clear locs
  locs=ceil(rand(num,2)*(N));
  return