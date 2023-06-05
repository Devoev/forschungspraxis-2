% Surface integration of normal component of a vectorfield. Here, a local
% coordinate system (u,v,w) is used. Such that this function works 
% properly,this system always needs to be a right-handed coordinate 
% system. Therefore, depending on the choice of the normal vector of the 
% surface, the local coordinates u and v have to be mapped to the 
% cartesian system accordingly.
% 
% msh       cartesian msh as created by cartMesh.
% vec       field vector 
% surf      description of integration surface
%           surf.ul/uh/vl/vh -> area of surface
%           ul/uh/vl/vh are indices (not coords!) in x-/y-/z-direction
%           surf.normal -> normal of surface (only positive entries!).
%           (does not work for surfaces that are not aligned with one of the axes)
%           The following coordinate mapping is induced:
%           surf.normal = [1,0,0] -> u := y, v := z
%			surf.normal = [0,1,0] -> u := x, v := z
%			surf.normal = [0,0,1] -> u := x, v := y
%           surf.w -> offset in normal direction
% 
% author: Thorben Casper
% created on: 2017/05/18

function val2D = intSurf(msh, vec, surf)

nx = msh.nx;
ny = msh.ny;
np = msh.np;
Mx = msh.Mx;
My = msh.My;
Mz = msh.Mz;

% selects the fit unit vector depending on the selected normal vector
Ndir = [[My;Mz;Mx],[Mx;Mz;My],[Mx;My;Mz]]*surf.normal';
% defines a scalar that takes care of the offset in w-direction and of
% chosing the correct component of the field vector vec
Ns = (surf.w-1)*Ndir(3)+[0,np,2*np]*surf.normal';

lengthUsurf = surf.uh-surf.ul+1;
lengthVsurf = surf.vh-surf.vl+1;

% initialize indices of primal facets (ipf) with offset in w-direction
ipf = Ns*ones(lengthVsurf,lengthUsurf);

% calculate u- and v-offsets as a matrix to be added to ipf
uOffsets = 0:Ndir(1):Ndir(1)*(lengthUsurf-1);
vOffsets = 0:Ndir(2):Ndir(2)*(lengthVsurf-1);
[Uoffset,Voffset] = meshgrid(uOffsets,vOffsets);

% add u- and v- offset to ipf
ipf = ipf + Uoffset + Voffset;
% flatten ipf and compute integral
ipf = ipf(:);
val2D = sum(vec(ipf));
