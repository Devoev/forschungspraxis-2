function [hnew,enew]=leapfrog(hold,eold,js,Mmui,Mepsi,c,dt)

% Berechnen der neuen magnetischen Spannung
 hnew = hold - dt*Mmui*c*eold;

% Berechnen der neuen elektrischen Spannung
 enew = eold + dt*Mepsi*(c' * hnew - js);

end

function [hnew,enew]=leapfrog(hold, eold, js, Mmui, Meps, c, Rmat, dt)

% Berechnen der neuen magnetischen Spannung
hnew = hold - dt*Mmui*c*eold;

% Berechnen der neuen elektrischen Spannung
enew = nullInv(nullInv(Rmat)+1/dt*Meps)*(1/dt*Meps*eold + c'*hnew - js);

end

function [hnew, enew] = leapfrog(h_old, e_old, e_given, h_given, Mmui, Mepsi, curl_matrix, dt, We, Wh, jbow)
    % LEAPFROG - Perform a time step using the leapfrog method for electromagnetic simulations.
    %
    %   [hnew, enew] = LEAPFROG(hold, eold, egiven, hgiven, Mmui, Mepsi, curl_matrix, dt, We, Wh)
    %
    %   Inputs:
    %     h_old        - Current magnetic field
    %     e_old        - Previous electric field
    %     e_given      - Boundary condition values for the electric field
    %                   (padded with zeros for domain DOFs)
    %     h_given      - Boundary condition values for the magnetic field
    %                   (padded with zeros for domain DOFs)
    %     Mmui        - Permeability inverse matrix
    %     Mepsi       - Permittivity matrix
    %     curl_matrix - Curl matrix for electromagnetic simulation
    %     dt          - Time step size
    %     We          - Electric field projector for boundary conditions
    %     Wh          - Magnetic field projector for boundary conditions
    %     jbow        - Given current densities
    %
    %   Outputs:
    %     hnew - New magnetic field after the time step
    %     enew - New electric field after the time step
    %
    %   Description:
    %   This function performs a time step using the leapfrog method for
    %   electromagnetic simulations. It updates the magnetic field 'hnew' and
    %   electric field 'enew' based on the input fields, material properties, and
    %   spatial operators. The leapfrog method is a time-stepping technique
    %   commonly used in computational electromagnetics.
    %
    %   We and Wh are projectors that enforce boundary conditions on the
    %   electric and magnetic fields, respectively. The 'egiven' and 'hgiven'
    %   inputs represent boundary condition values that are padded with zeros
    %   for degrees of freedom (DOFs) within the domain.
    %
    %   Example usage:
    %   [hnew, enew] = leapfrog(hold, eold, egiven, hgiven, Mmui, Mepsi, curl_matrix, dt, We, Wh);
    %
    %   See also: OTHER_RELATED_FUNCTIONS
    
    % Calculate the new magnetic field 'hnew'
    hnew = h_old - dt * Wh' * Mmui * curl_matrix * (e_given + We * e_old);

    % Calculate the new electric field 'enew'
    enew = e_old + dt * We' * Mepsi * (curl_matrix' * (Wh * hnew + h_given) - jbow);
end

