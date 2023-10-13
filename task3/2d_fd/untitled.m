A = [9 10 11 12;
     5 6 7 8;
     1 2 3 4;]
nx = 4;
ny = 3;
np = nx*ny;
all = 1:np; 

all(~(mod(all, nx)-3))