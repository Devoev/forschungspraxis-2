% one-dimensional integration of tangential component of vectorfield
% 
% Inputs:
% msh       cartesian msh as created by cartMesh.
% vec       field vector 
% line      description of integration path
%           line.u/v/w -> startpoint of integral
%           u/v/w are indices (not coordinates!) in x-/y-/z-direction
%           line.normal -> direction of integration path
%           e.g. [1,0,0], [0,1,0] or [0,0,1]
%           line.length -> number of points in integration path 
%
% author: Thorben Casper
% created on: 2017/05/18

function val1D = intEdge(msh, vec, line)

nx = msh.nx;   
ny = msh.ny;	
np = msh.np;
Mx=msh.Mx;
My=msh.My;
Mz=msh.Mz;

% find step size according to given normal vector (Mx,My or Mz)
Ndir = [Mx,My,Mz]*line.normal';
% find index of first primal edge to consider
ipeStart = line.u+(line.v-1)*My+(line.w-1)*Mz+[0,np,2*np]*line.normal';
% calculate the integral
val1D = sum(vec(ipeStart+(0:Ndir:(line.length-1)*Ndir)));
