//Airfoil
border upper(t=0, 1){x=t; y=0.17735*sqrt(t) - 0.075597*t - 0.212836*(t^2) + 0.17363*(t^3) - 0.06254*(t^4);}
border lower(t=1, 0){x = t; y=-(0.17735*sqrt(t) -0.075597*t - 0.212836*(t^2) + 0.17363*(t^3) - 0.06254*(t^4));}
border c(t=0, 2*pi){x=0.8*cos(t) + 0.5; y=0.8*sin(t);}
mesh Th = buildmesh(c(30) + upper(35) + lower(35));
plot(Th, ps="NACA0012.eps", bw=true);
