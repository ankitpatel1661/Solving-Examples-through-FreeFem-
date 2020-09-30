//Domain with Cardioid curve boundary
real b = 1, a = b;
border C(t=0, 2*pi){x=(a+b)*cos(t)-b*cos((a+b)*t/b); y=(a+b)*sin(t)-b*sin((a+b)*t/b);}
mesh Th = buildmesh(C(50));
plot(Th, ps="Cardioid.eps", bw=true);
