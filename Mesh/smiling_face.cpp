//Smiling face (Mouth is changeable)
real d=0.1; int m = 5; real a = 1.5, b = 2, c = 0.7, e = 0.01;

border F(t=0, 2*pi){x=a*cos(t); y=b*sin(t);}
border E1(t=0, 2*pi){x=0.2*cos(t)-0.5; y=0.2*sin(t)+0.5;}
border E2(t=0, 2*pi){x=0.2*cos(t)+0.5; y=0.2*sin(t)+0.5;}
func real st(real t){
    return sin(pi*t) - pi/2;
}
border C1(t=-0.5, 0.5){x=(1-d)*c*cos(st(t)); y=(1-d)*c*sin(st(t));}
border C2(t=0, 1){x=((1-d)+d*t)*c*cos(st(0.5)); y=((1-d)+d*t)*c*sin(st(0.5));}
border C3(t=0.5, -0.5){x=c*cos(st(t)); y=c*sin(st(t));}
border C4(t=0, 1){x=(1-d*t)*c*cos(st(-0.5)); y=(1-d*t)*c*sin(st(-0.5));}
border C0(t=0, 2*pi){x=0.1*cos(t); y=0.1*sin(t);}

mesh Th=buildmesh(F(10*m) + C1(2*m) + C2(3) + C3(2*m) + C4(3)
    + C0(m) + E1(-2*m) + E2(-2*m));
plot(Th, ps="SmileFace.eps", bw=true);
