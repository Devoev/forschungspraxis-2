% Path mesh functions
path_msh_func = '..\fit\mesh';
path_mat_func = '..\fit\matrices';
path_solver_func = '..\fit\solver';

% Add paths
addpath(path_msh_func, path_mat_func, path_solver_func)

elm = 100;
msh = cartMesh2D( linspace(1,elm,20*elm), linspace(1,elm,20*elm) );
np = msh.np;

[ds, dst, da, dat] = createGeoMats2D(msh);
[dsTE, dstTE, daTE, datTE] = createGeoMatsTE(msh);

normest(dsTE - ds(2*np+1:end, 2*np+1:end))
normest(dst(1:2*np, 1:2*np) - dstTE)
normest(da(1:2*np, 1:2*np) - daTE)
normest(dat(2*np+1:end, 2*np+1:end) - datTE)