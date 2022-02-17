load('HW6_G')

tsam = 1;
dsys = c2d(ss(G), tsam);
ad = dsys.A;
bd = dsys.B;
cd = dsys.C;