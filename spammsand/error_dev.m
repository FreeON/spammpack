<< NC`
<< NCAlgebra`


Clear[ a,b,c,d,f,g,h,i,j ];

SetNonCommutative[ Z,X,z,x ,e,dx]

z[0] := 1;
Z[0] := 1;
z[k_] := z[k-1]**m[ z[k-1]**s**z[k-1] ]
Z[k_] := (z[k-1]+tau*dz)**m[ Z[k-1]**s**Z[k-1] ]

D[ Z[3], tau]/.{tau->0}
