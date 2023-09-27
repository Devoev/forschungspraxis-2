
function [hnew,enew]=leapfrog_2D(hold,eold,js,Mmui,Mepsi,c,dt)
%% Description
%
% Leopfrog 2D time integration solver
%
% Input
% hold              h_bow last time step
% eold              e_bow last time step
% js                vector for source current (current time step)
% Mmui              inverse permeability matrix
% Mepsi             inverse permittivity matrix
% c                 primary curl matrix
% dt    	        value for time step
%
% Output
% hnew              h_bow for current time step
% enew              e_bow for current time step


%% Function definition

    % Berechnen der neuen magnetischen Spannung
    hnew = hold - dt*Mmui*c*eold;
    
    % Berechnen der neuen elektrischen Spannung
    enew = eold + dt*Mepsi*(c'*hnew - js);

end
