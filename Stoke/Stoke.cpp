/*********************************************************************
 * 
 * garanteed A posteroiri error estimate for Stokes problems using Uzawa algorithm
 * usage: FreeFem++ 
 * source article: https://who.rocq.inria.fr/Martin.Vohralik/Enseig/APost/a_posteriori.pdf  
 * page 111 for the Stokes problem
 *********************************************************************/

load "Element_Mixte" //For loading the RT1 element
system("clear");  


// define the mesh 
int numnode = 30;


macro Div(u1,u2) (dx(u1)+dy(u2)) //EOM
macro Grad(u) [dx(u),dy(u)]      //EOM


// mesh to define
mesh Th=square(numnode,numnode);
plot(Th,cmm="mesh",ps="mesh.eps");
real AreaDomain=Th.area;

// parameter in the Uzawa method, we define the value as 1 here
real alpha = 1; 

fespace Uh(Th,[P2,P2]);  // fonction space for  u=(u1,u2)
Uh [u1,u2],[v1,v2],[usave1,usave2];  
fespace Ph(Th,P1);       // function space for p
fespace Ph0(Th,P0);
Ph p, q;
Uh [uh1,uh2],[uvisu1,uvisu2];  
Ph ph,pvisu;

// Setting for the analytic solution
real beta = 0.44;    // to be change with each different case
func ue1=2*x^2*y*(x - 1)^2*(y - 1)^2 + x^2*y^2*(2*y - 2)*(x - 1)^2;
func ue2=- 2*x*y^2*(x - 1)^2*(y - 1)^2 - x^2*y^2*(2*x - 2)*(y - 1)^2;
func pe=x+y-1;
func dxue1=2*x*y^2*(2*y - 2)*(x - 1)^2 + 2*x^2*y*(2*x - 2)*(y - 1)^2 + x^2*y^2*(2*x - 2)*(2*y - 2) + 4*x*y*(x - 1)^2*(y - 1)^2;
func dyue1=2*x^2*y^2*(x - 1)^2 + 2*x^2*(x - 1)^2*(y - 1)^2 + 4*x^2*y*(2*y - 2)*(x - 1)^2;
func dxue2=- 2*x^2*y^2*(y - 1)^2 - 2*y^2*(x - 1)^2*(y - 1)^2 - 4*x*y^2*(2*x - 2)*(y - 1)^2;
func dyue2=- 2*x*y^2*(2*y - 2)*(x - 1)^2 - 2*x^2*y*(2*x - 2)*(y - 1)^2 - x^2*y^2*(2*x - 2)*(2*y - 2) - 4*x*y*(x - 1)^2*(y - 1)^2;
func f1=-24*x^4*y + 12*x^4 + 48*x^3*y - 24*x^3 - 48*x^2*y^3 + 72*x^2*y^2 - 48*x^2*y + 12*x^2 + 48*x*y^3 - 72*x*y^2 + 24*x*y - 8*y^3 + 12*y^2 - 4*y + 1;
func f2=48*x^3*y^2 - 48*x^3*y + 8*x^3 - 72*x^2*y^2 + 72*x^2*y - 12*x^2 + 24*x*y^4 - 48*x*y^3 + 48*x*y^2 - 24*x*y + 4*x - 12*y^4 + 24*y^3 - 12*y^2 + 1;


[uvisu1,uvisu2]=[ue1,ue2];
pvisu=pe;



// variational formulation to define:
varf a([u1,u2],[v1,v2])=int2d(Th)(Grad(u1)'*Grad(v1)+Grad(u2)'*Grad(v2))+int2d(Th)(f1*v1+f2*v2)
+on(1,2,3,4,u1=ue1,u2=ue2);
varf bc1([u1,u2],[v1,v2])=on(1,2,3,4,u1=1,u2=1);
varf b([u1,u2],[q])= int2d(Th)(Div(u1,u2)*q);
varf c([p],[q])=int2d(Th)(p*q);
matrix A=a(Uh,Uh,tgv=-1,solver=GMRES);
//set(A, solver=CG);
matrix B=b(Uh,Ph);
matrix Bt=B';
matrix C=c(Ph,Ph,solver=CG); // set solver=CG, to avoird the alert message  from C^-1; 


real[int] F=a(0,Uh,tgv=-1); //Boundary conditions
real[int] BC=bc1(0,Uh,tgv=-1); //Boundary conditions
Uh [b1,b2];


func real[int] DJ(real[int] &u)
{
real[int] Au=A*u;
return Au; // return of global variable ok
};


p=0; 
real epsgc= 1e-10;
// Uzawa iteration
int inumOuter=0;
while(inumOuter<10000){  //to be changed after	

  	b1[]=Bt*p[];  	b1[]+=F;   // second member for problem  Au=Btu+F  
	b1[]= BC ? F : b1[];
	u1[]= BC? F: u1[];
	
	LinearCG(DJ,u1[],b1[],veps=epsgc,nbiter=2000,verbosity=0);
	epsgc=-abs(epsgc); 

	real[int] Bu1=B*u1[];   
	Bu1*=alpha;	  
	real normBu1=sqrt(Bu1'*Bu1);
	if (normBu1<1e-10) break;		
		
	// update for p
	real[int] ptemp(p[].n);
	ptemp=C^-1*Bu1;  // C p^{k+1}=C p^{k} - alpha Bu
	p[]-=ptemp;
	p[]-=int2d(Th)(p)/AreaDomain;
	inumOuter++;
}
ph=p;
[uh1,uh2]=[u1,u2];     // save the solution of exact uzawa method
   
real errorgradu=sqrt(int2d(Th)( (Grad(u1)-[dxue1,dyue1])'*(Grad(u1)-[dxue1,dyue1])+(Grad(u2)-[dxue2,dyue2])'*(Grad(u2)-[dxue2,dyue2])  ));
real errorp=beta*sqrt(int2d(Th)( (p-pe)'*(p-pe) ));
real errortotal=errorgradu+errorp;
	

// Exact error computation	
fespace Whh(Th,P1);
varf indicatorErrorDiscu(unused,chiK) = int2d(Th)(chiK*((Grad(uh1)-[dxue1,dyue1])'*(Grad(uh1)-[dxue1,dyue1])+(Grad(uh2)-[dxue2,dyue2])'*(Grad(uh2)-[dxue2,dyue2]) ));
varf indicatorErrorDiscp(unused,chiK) = int2d(Th)(chiK*(beta*beta* (ph-pe)'*(ph-pe)   ));
Whh mapErrorDisc;
Whh mapErrorDiscu, mapErrorDiscp;
mapErrorDiscu[]=indicatorErrorDiscu(0,Whh);
mapErrorDiscp[]=indicatorErrorDiscp(0,Whh);
real GlobalErrorU=sqrt(mapErrorDiscu[].sum);
real GlobalErrorP=sqrt(mapErrorDiscp[].sum);
mapErrorDiscu=sqrt(mapErrorDiscu);
mapErrorDiscp=sqrt(mapErrorDiscp);
mapErrorDisc=mapErrorDiscu+mapErrorDiscp;


plot([uvisu1,uvisu2],value=1,fill=1,ps="analyticsolutionu.eps",cmm="Analytic solution u");
plot(pvisu,value=1,fill=1,ps="analyticsolutionp.eps",cmm="Analytic solution p");


plot([u1,u2],value=1,fill=1,ps="numericalsolutionu.eps",cmm="Numerical solution uh");
plot(p,value=1,fill=1,ps="numericalsolutionp.eps",cmm="Numerical solution ph");


plot(mapErrorDiscu,ps="ErrorDiscU.eps",value=1, fill=1,cmm="ErrorDiscretizationU");	
plot(mapErrorDiscp,ps="ErrorDiscP.eps",value=1, fill=1,cmm="ErrorDiscretizationP");	
plot(mapErrorDisc,ps="ErrorDiscTotal.eps",value=1, fill=1,cmm="ErrorDiscretization");	


 
//  A posteriori estimator
mesh ThaI=Th;
fespace Pa(ThaI,P1dc);
fespace Xh(Th,P1);
Xh phia=0;
mesh[int] Tha(Th.nv);
// construction of the patch for each node on the mesh
for(int i=0; i <Th.nv; ++i)
{
	phia[][i]=1; // phia 
	Tha[i] = trunc(Th,phia>0,label=10);
	if (Th(i).label==0) Tha[i]=change(Tha[i],flabel=10); 
	phia[][i]=0; // phia 
}

fespace VPa(ThaI,[RT1,RT1,P1dc,P1dc]); 
fespace Vh(Th,[RT1,RT1]);
real[int] areaPatcha(Th.nv);
matrix[int] Aa(Th.nv);
real[int][int] Fa(Th.nv);
Vh [d11,d12,d21,d22]=[0,0,0,0];
Vh [d11p,d12p,d21p,d22p]=[0,0,0,0];
real epsregua = numnode*numnode*1e-10;   

for( int i=0; i < Th.nv ;i++) 
{
	phia[][i]=1; // phia 
	ThaI = Tha[i];
	Pa Chia =1; 
	VPa [d11t,d12t,d21t,d22t,qa1,qa2],[vh11,vh12,vh21,vh22,phih1,phih2];
	varf Varflocala([d11t,d12t,d21t,d22t,qa1,qa2],[vh11,vh12,vh21,vh22,phih1,phih2]) = 
	int2d(ThaI)( [d11t,d12t,d21t,d22t]'*[vh11,vh12,vh21,vh22]
	+ (qa1*Div(vh11,vh12)+qa2*Div(vh21,vh22))
	+ (Div(d11t,d12t)*phih1+Div(d21t,d22t)*phih2)
	- (epsregua*(qa1*phih1+qa2*phih2))) // regularization;
	+ on(10,d11t=0,d12t=0,d21t=0,d22t=0);
	matrix Art=Varflocala(VPa,VPa);	
	Aa[i]= Art; 
	set(Aa[i],solver=sparsesolver);

	phia[][i]=0;
}

int[int][int] I(Th.nv); // number for Vh to Vhi
real[int][int] epsI(Th.nv); 
Vh [num,num1,num2,num3];
num[]=1:num[].n; 
for(int i=0; i <Th.nv; ++i)
{ 
   ThaI = Tha[i];
   VPa [ul1,ul2,ul3,ul4,pl1,pl2];
   [ul1,ul2,ul3,ul4,pl1,pl2] = [num,num1,num2,num3,1e9,1e9];
   I[i].resize(VPa.ndof); 
   epsI[i].resize(VPa.ndof); 
   for(int j=0;j<ul1[].n;++j) {
	  epsI[i][j]= ul1[][j] < 0 ? -1 : 1;
	  I[i][j]= abs(ul1[][j])-0.5 ;  
	  if(I[i][j]>= 999999999)I[i][j]=-1;
	}
}
Vh [d11res,d12res,d21res,d22res]=[0,0,0,0];

for(int i=0; i < Th.nv ;++i) 
{
	phia[][i]=1; 
	ThaI = Tha[i];	
     	Pa Chia =1; 
	fespace ProRT1(ThaI,[RT1,RT1]);
	ProRT1 [pd11,pd12,pd21,pd22];
	[pd11,pd12,pd21,pd22]=([dx(u1),dy(u1),dx(u2),dy(u2)]-[p,0,0,p])*phia;
	
	VPa [d11tmp,d12tmp,d21tmp,d22tmp,qa1,qa2],[v11,v12,v21,v22,phih1,phih2];	
       
	varf l([d11t,d12t,d21t,d22t,qa1,qa2],[vh11,vh12,vh21,vh22,phih1,phih2]) 
	= int2d(ThaI)([pd11,pd12,pd21,pd22] '*[vh11,vh12,vh21,vh22]
	  -(f1*phia*phih1+f2*phia*phih2)
          -(p*phih1*dx(phia)+p*phih2*dy(phia))
	  +(Grad(u1)'*Grad(phia)*phih1+Grad(u2)'*Grad(phia)*phih2))
	  + on(10,d11t = 0, d12t = 0, d21t=0, d22t=0);

	VPa [F,F1,F2,F3,F4,F5] ; 
	F[] = l(0,VPa);
	Fa[i].resize(F[].n);
	Fa[i]=F[];
	d11tmp[] = Aa[i]^-1*Fa[i];		
	real[int] so(d11tmp[].n);
	so=d11res[](I[i]);
	d11tmp[]= d11tmp[].*epsI[i];
	d11tmp[] += so;
	d11res[](I[i]) = d11tmp[];		
	phia[][i]=0;	
}	
[d11,d12,d21,d22]=[d11res,d12res,d21res,d22res];


Ph0 fh1,fh2;
fh1=f1;
fh2=f2;
// computation for the estimator
varf indicatorEtaDisc(unused,chiK) = int2d(Th)(chiK*(Div(u1,u2)'*Div(u1,u2))/beta/beta)  
+  int2d(Th)(chiK*hTriangle*hTriangle/pi/pi*((f1-fh1)*(f1-fh1)+(f2-fh2)*(f2-fh2))) 
+  int2d(Th)(chiK*(  ([dx(u1),dy(u1),dx(u2),dy(u2)]-[d11,d12,d21,d22]-[p,0,0,p])'*([dx(u1),dy(u1),dx(u2),dy(u2)]-[d11,d12,d21,d22]-[p,0,0,p]))); 	     	     
Whh mapEtaDisc;
mapEtaDisc[]=indicatorEtaDisc(0,Whh);
real GlobalEstimator=sqrt(mapEtaDisc[].sum);
mapEtaDisc=sqrt(mapEtaDisc);
plot(mapEtaDisc,  ps="MapDiscretisationEstimator.eps",value=1, fill=1,cmm="DiscretizationEstimator");

cout<<"************************************ "<<endl;
cout<<"GlobalError on u= " <<GlobalErrorU<<endl;
cout<<"GlobalError on p= " <<GlobalErrorP<<endl;
cout<<"GlobalError = " <<GlobalErrorU+GlobalErrorP<<endl;
cout<<"GlobalEstimator = " <<GlobalEstimator<<endl;
