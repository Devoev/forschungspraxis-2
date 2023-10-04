function [MAT, bc] = generate_MAT(msh, bc, material_regions, operators)
% generate_MAT Generates all calculation relevant matrices -> Operators and
% material matrices and store them in one object MAT
%
% Inputs:
%   msh                 - Mesh struct
%   bc                  - Boundary condition struct
%   material_regions    - Struct containing information about the material
%                         regions
%                         Example for regions for relative permittivity:
%                         boxesEpsilonR(1).box = [1, nx, 1, ny];
%                         boxesEpsilonR(1).value = 1;
%                         material_regions.boxesEpsilonR = boxesEpsilonR;
%                         -> Same ist needed for inverse permeability
%                         
%   operators           - List of desired operators like
%                         operators = ["CurlP", "SourceD"];
%                         -> Only operators named in the list are returned
%                         in MAT object. 
%                         CurlP -> Primary curl
%                         SourceP -> Primary source
%                         SourceD -> Dual source
%
% Outputs:
%   MAT                 - Struct containing all edited matrices


%% Add basic constants to boundary condition object bc and MAT

bc.epsilon0 = material_regions.epsilon0;
bc.mu0i = material_regions.mu0i;

MAT.epsilon0 = material_regions.epsilon0;
MAT.mu0i = material_regions.mu0i;


%% Generate matrices for calculation

% Create curl, source and geometric matirces
[c, s, st] = createTopMats_2D(msh);
[ds, dst, da, dat] = createGeoMats_2D(msh);

% Create permittivity matrix if associated regions are defines
if any(ismember(fieldnames(material_regions),"boxesEpsilonR"))
    rel_eps = boxMesher_2D(msh, material_regions.boxesEpsilonR, material_regions.epsilon0);
    meps = createMeps_2D(msh, ds, da, dat, rel_eps, material_regions.epsilon0);
    MAT.meps = meps;
end

% Create inverse permeability matrix if associated regions are defines
if any(ismember(fieldnames(material_regions),"boxesMuiR"))
    rel_mui = boxMesher_2D(msh, material_regions.boxesMuiR, material_regions.mu0i);
    mmui = createMmui_2D(msh, ds, dst, da, rel_mui, material_regions.mu0i);
    MAT.mmui = mmui;
end

% Add desired operator matrices to output matrix object MAT
if any(ismember(operators,"CurlP"))
    MAT.c = c;
end
if any(ismember(operators,"SourceP"))
    MAT.s = s;
end
if any(ismember(operators,"SourceD"))
    MAT.st = st;
end


end 
