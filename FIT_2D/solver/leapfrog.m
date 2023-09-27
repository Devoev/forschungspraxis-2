function [hnew,enew]=leapfrog(hold,eold,js,Mmui,Mepsi,c,dt)

% Berechnen der neuen magnetischen Spannung
hnew = hold - dt*Mmui*c*eold;

% Berechnen der neuen elektrischen Spannung
enew = eold + dt*Mepsi*(c'*hnew - js);

end
