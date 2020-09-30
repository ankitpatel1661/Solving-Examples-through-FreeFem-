/*********************************************************************
 * 
 * garanteed A posteroiri error estimate for laplace problem
  
 *********************************************************************/

load "Element_Mixte" //For loading the RT1 element
system("clear");  

   
// Table for computational time
real[int] TimeSave(20);
int timeConsPatchThaI=1;
int timeConsGlobMat=2;
int timeCompEtaAlg=3;
int timeMapEtaAlg=4; 
int timeCompEtaDisc=5;
int timeMapEtaDisc=6;
int timeConsRHS=7;
int timeSolvRHS=8;
int timeCompEta=9;
int timeMapEta=10;
int timeSover=11;
int timeProRHS=12;
int timeCompEtaPart1=13;
int timeCompEtaPart2=14;
int timeCompEtaPart3=15;
int timeCompEtaPart4=16;
int timeCompEtaPart5=17;
int timeProEta=18;
int timeCompEtaOptmised=19;
 

// define the mesh 
int numnode = 80;
mesh Th = square(numnode,numnode);



// case with an analytic solution, and its derivatives [dx(u),dy(u)]
func uExact = x*(x-1)*y*(y-1);    
func dxuExact=(2*x-1)*y*(y-1);
func dyuExact=(2*y-1)*x*(x-1); 
// boundary condition 
func gd=x*(x-1)*y*(y-1);   
func f = -2*(x*x+y*y)+2*(x+y);

// macro definition
macro Div(u1,u2) (dx(u1)+dy(u2)) //EOM
macro Grad(u) [dx(u),dy(u)]      //EOM


// functional space for solution 
fespace Vh(Th,P1);
Vh uh, vh;
Vh uvisu=uExact;   //visualisation of the exact solution, projection in Vh space


// variational formulation
varf a(uh,vh)=int2d(Th)(Grad(uh)'*Grad(vh)) + int2d(Th)(f*vh) + on(1,2,3,4,uh=gd);

// solve the laplace problem exactly
matrix A=a(Vh,Vh,solver=CG);
real[int] b=a(0,Vh);
real TimeCurrent=clock();
uh[]=A^-1*b;
TimeSave[timeSover]=clock()-TimeCurrent;





//visualisation of the analytic solution and numerical solution
plot(uvisu,value=1,fill=1,ps="analyticsolution.eps",cmm="analytic solution");
plot(uh,value=1,fill=1,ps="numericalsolution.eps",cmm="numerical solution");



real GlobalError;
func int ErrorMapGeneration()
{
    cout<<" TO generate the distribution of error map ... "<<endl;
    // computation for the exact error and map of error distribution
    fespace Wh(Th,P0); // constant in each element
    varf indicatorErrorTotal(unused,chiK) = int2d(Th)(chiK*((Grad(uh)-[dxuExact,dyuExact])'*(Grad(uh)-[dxuExact,dyuExact])));
    Wh mapErrorTotal;
    mapErrorTotal[]=indicatorErrorTotal(0,Wh);
    GlobalError=sqrt(mapErrorTotal[].sum);
    mapErrorTotal=sqrt(mapErrorTotal);
    plot(mapErrorTotal, ps="MapDiscretisationError.eps",value=1, fill=1,cmm="DiscretisationError");
}



//************************ Computational part: a posteriori estimator  **************************
func int versionclassic()
{
mesh ThaI=Th;
fespace Pa(ThaI,P1dc);
fespace VPa(ThaI,[RT1,P1dc]); 
fespace VRT(Th,RT1);


Vh phia=0;
real epsregua = numnode*numnode*1e-10;   //regularization for sub-problem,
mesh[int] Tha(Th.nv);
matrix[int] Aa(Th.nv);


// construction of the patch for each node on the mesh
TimeCurrent=clock();
for(int i=0; i <Th.nv; ++i)
{
    phia[][i]=1;
    Tha[i] = trunc(Th,phia>0,label=10);
    if (Th(i).label==0) Tha[i]=change(Tha[i],flabel=10); 
    phia[][i]=0;
}
TimeSave[timeConsPatchThaI]=clock()-TimeCurrent;


// construction of the local problem
TimeCurrent=clock();
for( int i=0; i < Th.nv ;++i) 
{
//  to bluild the basic function 
    phia[][i]=1; // phia 
    ThaI = Tha[i];
    Pa Chia =1; 
    VPa [s1,s2,ra],[v1,v2,q];
    //matrix	
    varf a([s1,s2,ra],[v1,v2,q]) = int2d(ThaI)( [s1,s2]'*[v1,v2] - ra*Div(v1,v2) - Div(s1,s2)*q - epsregua*ra*q )+ on(10, s1=0,s2=0);
    matrix Art=a(VPa,VPa);
    Aa[i]= Art; 		
    set(Aa[i],solver=sparsesolver);
    phia[][i]=0;
}
TimeSave[timeConsGlobMat]=clock()-TimeCurrent;


// script de F. Hecht:  number for global matrix to local matrix
int[int][int] I(Th.nv);  
real[int][int] epsI(Th.nv); 
VRT [num,num1];
real[int][int] Fa(Th.nv);
num[]= 1: num[].n ; //  
for(int i=0; i <Th.nv; ++i)
{ 
    ThaI = Tha[i];
    VPa [ul,vl,p];
    [ul,vl,p] = [num,num1,1e9];
    I[i].resize(VPa.ndof); 
    epsI[i].resize(VPa.ndof); 
    for(int j=0;j<ul[].n;++j) 
    {
	epsI[i][j]= ul[][j] < 0 ? -1 : 1;
	I[i][j]= abs(ul[][j])-0.5 ;  
	if(I[i][j]>= 999999999)I[i][j]=-1;
    } 
}


//for loop computing right hand side and solve local problem
VRT [sigma1,sigma2]=[0,0];
TimeCurrent=clock();
for(int i=0; i < Th.nv ;++i) 
{
    phia[][i]=1; // phia 
    ThaI = Tha[i];
    Pa Chia =1; 
    VPa [s1,s2,ra],[v1,v2,q];
    //right hand side			
    varf l([uv1,uv2,uq],[v1,v2,q]) = int2d(ThaI)( - phia*(Grad(uh)'*[v1,v2] + f*q) + Grad(phia)'*Grad(uh)*q)
				    + on(10,uv1 = 0, uv2 = 0);
    VPa [F,F1,F2] ; F[] = l(0,VPa);
    Fa[i].resize(F[].n);
    Fa[i]=F[]; 
    s1[] = Aa[i]^-1*Fa[i];
    real[int] so(s1[].n);
    so=sigma1[](I[i]);
    s1[]= s1[].*epsI[i];
    s1[] += so;
    sigma1[](I[i]) = s1[];
    phia[][i]=0;
}
TimeSave[timeConsRHS]=clock()-TimeCurrent;


real GlobalEstimator;
func int EstimatorMapGeneration()
{
    TimeCurrent=clock();
    cout<<" TO generate the distribution of estimator map ... "<<endl;
    // computation for the estimator
    fespace Wh(Th,P0); // constant in each element
    varf indicatorEtaDisc(unused,chiK) = int2d(Th)(chiK*( (Grad(uh) + [sigma1,sigma2]) '* (Grad(uh) + [sigma1,sigma2])  ))
				+ int2d(Th)(chiK*(hTriangle*hTriangle/pi/pi* (f-Div(sigma1,sigma2))*(f-Div(sigma1,sigma2)) ));
    Wh mapEtaDisc;
    mapEtaDisc[]=indicatorEtaDisc(0,Wh);
    GlobalEstimator=sqrt(mapEtaDisc[].sum);
    mapEtaDisc=sqrt(mapEtaDisc);
    //mapEtaDisc=mapEtaDisc/mapEtaDisc[].max;
    plot(mapEtaDisc,  ps="MapDiscretisationEstimator.eps",value=1, fill=1,cmm="DiscretizationEstimator");
    TimeSave[timeMapEta]=clock()-TimeCurrent;
}

TimeCurrent=clock();
// computation for the estimator
real resEstimator=int2d(Th)(( (Grad(uh) + [sigma1,sigma2]) '* (Grad(uh) + [sigma1,sigma2])  ))
	     + int2d(Th)((hTriangle*hTriangle/pi/pi* (f-Div(sigma1,sigma2))*(f-Div(sigma1,sigma2)) ));
resEstimator=sqrt(resEstimator);
//mapEtaDisc=mapEtaDisc/mapEtaDisc[].max;
TimeSave[timeCompEta]=clock()-TimeCurrent;

ErrorMapGeneration();
EstimatorMapGeneration();
cout<<"GlobalError = " <<GlobalError<<endl;
cout<<"GlobalEstimator = " <<GlobalEstimator<<endl;  


}
 

real TimeVersionClassic=clock();
versionclassic();
TimeVersionClassic=clock()-TimeVersionClassic;

cout<<"Time for computation estimator  (without map)= " <<TimeSave[timeCompEta]<< endl;
cout<<"Time for computation estimator (with map) = " <<TimeSave[timeMapEta]<< endl;
cout<<"Time for right-hand-side of local problem = " <<TimeSave[timeConsRHS]<< endl;

cout<<" ***************************" <<endl;
cout<<"Total time for solving Ax=b = " <<TimeSave[timeSover]<<endl;
cout<<"Time for constructions of all the patches (Pre Proceeding)  = " <<TimeSave[timeConsPatchThaI]<<endl;
cout<<"Time for constructions of the global matrix (Pre Proceeding)= "<<TimeSave[timeConsGlobMat]<< endl;



 






