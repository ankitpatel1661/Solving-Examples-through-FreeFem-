//Cool air (green) comes from the lower left and mix with hot air (magenta), the right boundary is free. This is Navier-Stokes-Boussinesq integrated with P1-bubble P1 mixte finite element.
//A very small example 2d of how to solve the Poisson equation on a L shape

border aaa(t=0,1){x=t;y=0;};
border bbb(t=0,0.5){x=1;y=t;};
border ccc(t=0,0.5){x=1-t;y=0.5;};
border ddd(t=0.5,1){x=0.5;y=t;};
border eee(t=0.5,1){x=1-t;y=1;};
border fff(t=0,1){x=0;y=1-t;};
mesh Th = buildmesh (aaa(6) + bbb(4) + ccc(4) +ddd(4) + eee(4) + fff(6));
fespace Vh(Th,P1);   //  to change P1 in P2 to make P2 finite element.
Vh u=0,v;
func f= 1;
func g= 0;
int i=0;
real error=0.1, coef= 0.1^(1./5.);
problem Probem1(u,v,solver=CG,eps=-1.0e-6) =
    int2d(Th)(  dx(u)*dx(v) + dy(u)*dy(v)) 
  + int2d(Th) ( v*f ) 
  + on(aaa,bbb,ccc,ddd,eee,fff,u=g)  ;
  
for (i=0;i< 10;i++)
{   
  real d = clock();
  Probem1; //  solve the problem 
  plot(u,Th,wait=1);
  Th=adaptmesh(Th,u,inquire=1,err=error);
  error = error * coef;
} ;
