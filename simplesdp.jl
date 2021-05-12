


tmax = 100;
xmax = 10;


f = 0.2;
p = 1;

S = zeros(xmax,tmax);

S[2:xmax,tmax] .= 1;

for i=collect((tmax-1):-1:1)
    for x = xc:xmax
        