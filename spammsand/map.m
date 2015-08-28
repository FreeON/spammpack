f[x_]:= x*scalshift+shftshift
g[x_]:= x*scalmapp +shftmapp

  h=Simplify[Expand[g[f[x]] ]];

Print[D[h,x]]

