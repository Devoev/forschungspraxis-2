function [ebow_new,hbow_new]=solve_FullLeapfrog_2d_td(ebow_old,hbow_old,e_exi_old,e_exi_new,jsbow,mmui,mepsi,kaps,c,dt,W)
% solve_FullLeapfrog_2d_td executes an explicit time integration step 
% utilizing the leapfrog algorithmn including a conductivity matrix
%
% Inputs:
%   ebow_old    - Vector with integrated E-field values of start time 
%                 size (3*np x 1)
%   hbow_old    - Vector with integrated H-field values of start time
%                 size (3*np x 1)
%   e_exi_old   - Vector with integrated E-field excitation values of 
%                 start time size (3*np x 1)
%   e_exi_new   - Vector with integrated E-field excitation values of 
%                 end time t+dt size (3*np x 1)
%   jsbow       - Vector with integrated current excitation values of 
%                 start time size (3*np x 1)
%   mmui        - Inverse permeability matrix
%   mepsi       - Inverse permittivity matrix
%   kaps        - Conductivity matrix
%   c           - Primary curl matrix
%   dt          - Size of timestep (double value)
%   W           - Projector matrix W
%
% Outputs:
%   ebow_new    - vector conatining new grid voltages (3*np x 1)
%   hbow_new    - vector conatining new integrated H-field values 
%                 (3*np x 1)


    % Berechnen der neuen magnetischen Spannung
    hbow_new = hbow_old - dt*mmui*c* (e_exi_old + W*W'*ebow_old);
    
    % Berechnen der neuen elektrischen Spannung
    ebow_new = (e_exi_old + W*W'*ebow_old) + dt*mepsi*(c'*hbow_new - kaps * ebow_old - jsbow);
    ebow_new = e_exi_new + W*W'*ebow_new;


end
