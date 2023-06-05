function [hnew,enew]=leapfrog(hold, eold, js, Mmui, Meps, c, Rmat, dt)

% Berechnen der neuen magnetischen Spannung
hnew = hold - dt*Mmui*c*eold;

% Berechnen der neuen elektrischen Spannung
enew = nullInv(nullInv(Rmat)+1/dt*Meps)*(1/dt*Meps*eold + c'*hnew - js);

end