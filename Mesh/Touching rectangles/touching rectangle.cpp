//Two touchinging rectangles

border a(t=0, 1){x=t; y=0;};
border b(t=0, 1){x=1; y=t;};
border c(t=1, 0){x=t; y=1;};
border d(t=1, 0){x=0; y=t;};
border c1(t=0, 1){x=t; y=1;};
border e(t=0, 0.2){x=1; y=1+t;};
border f(t=1, 0){x=t; y=1.2;};
border g(t=0.2, 0){x=0; y=1+t;};
int n=1;
mesh th = buildmesh(a(10*n) + b(10*n) + c(10*n) + d(10*n));
mesh TH = buildmesh(c1(10*n) + e(5*n) + f(10*n) + g(5*n));
plot(th, TH, ps="TouchSide.eps");
