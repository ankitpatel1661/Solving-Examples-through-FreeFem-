//Mesh without and with hole 
border a(t=0, 2*pi){x=cos(t); y=sin(t); label=1;}
border b(t=0, 2*pi){x=0.3+0.3*cos(t); y=0.3*sin(t); label=2;}
plot(a(50) + b(30)); //to see a plot of the border mesh
mesh Thwithouthole = buildmesh(a(50) + b(30));
mesh Thwithhole = buildmesh(a(50) + b(-30));
plot(Thwithouthole, ps="Thwithouthole.eps")
plot(Thwithhole, ps="Thwithhole.eps")
