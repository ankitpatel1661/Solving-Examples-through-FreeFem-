// Mesh
mesh Th = square(3, 3);

// Fespace
fespace Vh(Th, P1);
Vh u=0;

// Loop on all degrees of freedom
int n=u.n;
for (int i = 0; i < n; i++){
    u[][i] = 1; // The basis function i
    plot(u, wait=true);
    mesh Sh1 = trunc(Th, abs(u)>1.e-10, split=5, label=2);
    plot(Th, Sh1, wait=true, ps="trunc"+i+".eps");
    u[][i] = 0; // reset
}
