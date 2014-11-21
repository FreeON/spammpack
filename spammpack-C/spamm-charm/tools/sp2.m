function P = sp2(F)
  bounds = gershgorin(F);
  P = (eye(size(F,1))*bounds(2))/(bounds(2)-bounds(1));
  for i = 1, 100
    P2 = P*P;
    traceP2 = trace(P2);
    if(trace(P2)-norm < 2*trace(P)-trace(P2)-norm)
      P = P2;
    else
      P = 2*P-P2;
    end
