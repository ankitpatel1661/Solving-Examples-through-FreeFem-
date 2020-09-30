//Domain with Cassini egg curve boundary
border C(t=0, 2*pi) {x=(2*cos(2*t)+3)*cos(t); y=(2*cos(2*t)+3)*sin(t);}
mesh Th = buildmesh(C(50));
plot(Th, ps="Cassini.eps", bw=true);

//U-shaped
real d = 0.1; //width of U-shape
border L1(t=0, 1-d){x=-1; y=-d-t;}
border L2(t=0, 1-d){x=-1; y=1-t;}
border B(t=0, 2){x=-1+t; y=-1;}
border C1(t=0, 1){x=t-1; y=d;}
border C2(t=0, 2*d){x=0; y=d-t;}
border C3(t=0, 1){x=-t; y=-d;}
border R(t=0, 2){x=1; y=-1+t;}
border T(t=0, 2){x=1-t; y=1;}
int n = 5;
mesh Th = buildmesh(L1(n/2) + L2(n/2) + B(n) + C1(n) + C2(3) + C3(n) + R(n) + T(n));
plot(Th, ps="U-shape.eps", bw=true);

//V-shape cut
real dAg = 0.02; //angle of V-shape
border C(t=dAg, 2*pi-dAg){x=cos(t); y=sin(t);};
real[int] pa(2), pb(2), pc(2);
pa[0] = cos(dAg);
pa[1] = sin(dAg);
pb[0] = cos(2*pi-dAg);
pb[1] = sin(2*pi-dAg);
pc[0] = 0;
pc[1] = 0;
border seg1(t=0, 1){x=(1-t)*pb[0]+t*pc[0]; y=(1-t)*pb[1]+t*pc[1];};
border seg2(t=0, 1){x=(1-t)*pc[0]+t*pa[0]; y=(1-t)*pc[1]+t*pa[1];};
mesh Th = buildmesh(seg1(20) + C(40) + seg2(20));
plot(Th, ps="V-shape.eps", bw=true);
