function plotEdgeVoltage(msh,volt,izcut,bc)
% plotEdgeVoltage(msh,volt,izcut,bc) plots the field corresponding to edge voltage vector volt with the help of quiver and quiver3
% 
% msh   = struct with grid data; must contain at least:
%         msh.nx, msh.ny, msh.nz
%         msh.xmesh, msh.ymesh, msh.zmesh
% volt  = edge voltage (not the field itself!)
% izcut = index of zcut-plane (2D-mode)
%         or array of indices (e.g. [1:nz]) (3D-mode)
%         (optional, default: 1:nz)
% bc    = array of boundary conditions: [xlow xhigh ylow yhigh zlow zhigh]
%         1 = electric
%         0 = continuous
%		  (optional, default: [1 1 1 1 1 1])
%
% author: unknown
% editor: Thorben Casper

if nargin<4
    bc = [1 1 1 1 1 1];
end
if nargin<3
   izcut = 1:msh.nz;
end

ratio = 1;

multiple = numel(izcut) > 1;
if multiple
    narrows = 300;
    narr2d = round(narrows/numel(izcut));
else
    narr2d = 200;
end
xm = reshape(msh.xmesh, msh.nx,1);
ym = reshape(msh.ymesh, msh.ny,1);
zm = reshape(msh.zmesh, msh.nz,1);

x1 = xm(1); x2 = xm(end);
y1 = ym(1); y2 = ym(end);
z1 = zm(1); z2 = zm(end);

xx = (x2-x1);
yy = (y2-y1);
zz = (z2-z1);
if ~ratio
%    stretch = [1 (sqrt(5)+1)/2 *yy/xx];
    if zz~=0
        stretch = [1 xx/yy xx/zz];
    else
        stretch = [1 xx/yy 1];
    end
    
    xm = stretch(1)*xm;
    ym = stretch(2)*ym;
    zm = stretch(3)*zm;
else
    stretch = [1 1 1];
end

nxplot = round(sqrt(narr2d*(stretch(1)*xx)/(stretch(2)*yy)));
nyplot = round(narr2d/nxplot);

nx = msh.nx; ny = msh.ny; nz = msh.nz;

[~,~,izgrid] = meshgrid(1:nx, 1:ny, 1:nz);
xgrid = repmat(reshape(xm,[1,nx,1]),[ny,1,nz]);
ygrid = repmat(reshape(ym,[ny,1,1]),[1,nx,nz]);
zgrid = repmat(reshape(zm,[1,1,nz]),[ny,nx,1]);

% shift grid (allocation of components)
xdual = (xm + [xm(2:nx); xm(end)])/2;
ydual = (ym + [ym(2:ny); ym(end)])/2;
zdual = (zm + [zm(2:nz); zm(end)])/2;
xgrid_x = repmat(reshape(xdual,[1,nx,1]),[ny,1,nz]);
ygrid_y = repmat(reshape(ydual,[ny,1,1]),[1,nx,nz]);
zgrid_z = repmat(reshape(zdual,[1,1,nz]),[ny,nx,1]);

vec3 = reshape(volt, length(volt),1);

ddx = [reshape(diff(msh.xmesh),msh.nx-1,1); 0];
ddy = [reshape(diff(msh.ymesh),msh.ny-1,1); 0];
ddz = [reshape(diff(msh.zmesh),msh.nz-1,1); 0];
mat_dsx = repmat(ddx,msh.ny*msh.nz,1);
mat_dsy = repmat(reshape(repmat(ddy',msh.nx,1), msh.nx*msh.ny, 1),msh.nz,1);
mat_dsz = reshape(repmat(ddz,1,msh.nx*msh.ny)', msh.nx*msh.ny*msh.nz, 1);

fff = find(mat_dsx);
mat_dsx_i = mat_dsx; mat_dsx_i(fff) = 1./mat_dsx(fff);
fff = find(mat_dsy);
mat_dsy_i = mat_dsy; mat_dsy_i(fff) = 1./mat_dsy(fff);
fff = find(mat_dsz);
mat_dsz_i = mat_dsz; mat_dsz_i(fff) = 1./mat_dsz(fff);

if length(vec3) == 2*nx*ny
    vec3 = [vec3; zeros(nx*ny,1)];  % enlarge 2D to 3D vector
end

vec3 = vec3 .* [mat_dsx_i; mat_dsy_i; mat_dsz_i];

vx = reshape(vec3(           1:  nx*ny*nz), [nx,ny,nz]);
vy = reshape(vec3(  nx*ny*nz+1:2*nx*ny*nz), [nx,ny,nz]);
vz = reshape(vec3(2*nx*ny*nz+1:3*nx*ny*nz), [nx,ny,nz]);

vlen = size(vx,1)*size(vx,2)*numel(izcut);
vxmax = max(reshape(abs(vx(:,:,izcut)),[vlen,1]));
vymax = max(reshape(abs(vy(:,:,izcut)),[vlen,1]));
vzmax = max(reshape(abs(vz(:,:,izcut)),[vlen,1]));
vmax = max([vxmax vymax vzmax]);

% plotgrid 1/2 Schritt vom Rand weg: also /nxplot statt /(nxplot-1)
dxplot = (xm(end)-xm(1))/nxplot;
dyplot = (ym(end)-ym(1))/nyplot;
[ixplot,iyplot] = meshgrid(1:nxplot, 1:nyplot);
xplot = xm(1) + (ixplot-1/2)*dxplot;
yplot = ym(1) + (iyplot-1/2)*dyplot;

% scale = min(dxplot,dyplot)/vmax;
scale = 0.5*min(dxplot,dyplot)/vmax;

hhh = ishold;
if ~hhh
    clf
end
hold on

for iz = izcut
    zplot = zm(iz) * ones(size(xplot));

    % in meshgrid wird x und y andersrum angeordnet, muss deshalb hier getauscht werden
    vxcut = vx(:,:,iz)';
    vycut = vy(:,:,iz)';
    if iz==1
        vzcut = vz(:,:,iz)';
    else
        vzcut = (vz(:,:,iz)' + vz(:,:,iz-1)')/2;
    end

    % extend grid at lower boundary for proper interpolation
    xgrid_x_cut = [xm(1)*ones(ny,1) xgrid_x(:,:,iz)];
    xgrid_x_cut(:,nx+1) = xm(end)*ones(ny,1);
    ygrid_y_cut = [ym(1)*ones(1,nx); ygrid_y(:,:,iz)];
    ygrid_y_cut(ny+1,:) = ym(end)*ones(1,nx);

    % new boundary components: copy from next slice or set to zero
    if bc(1)
        vxcut = [vxcut(:,1) vxcut];
    else
        vxcut = [zeros(ny,1) vxcut];
    end
    if bc(2)
        vxcut(:,nx+1) = vxcut(:,nx);
    else
        vxcut(:,nx+1) = zeros(ny,1);
    end
    if bc(3)
        vycut = [vycut(1,:); vycut];
    else
        vycut = [zeros(1,nx); vycut];
    end
    if bc(4)
        vycut(ny+1,:) = vycut(ny,:);
    else
        vycut(ny+1,:) = zeros(1,nx);
    end

    ygrid_x_cut = ygrid(:,:,iz); ygrid_x_cut = [ygrid_x_cut(:,1) ygrid_x_cut];
    xgrid_y_cut = xgrid(:,:,iz); xgrid_y_cut = [xgrid_y_cut(1,:); xgrid_y_cut];

    wx = scale*griddata(xgrid_x_cut  , ygrid_x_cut    , vxcut ,xplot,yplot);
    wy = scale*griddata(xgrid_y_cut  , ygrid_y_cut    , vycut ,xplot,yplot);
    wz = scale*griddata(xgrid(:,:,iz), ygrid(:,:,iz)  , vzcut ,xplot,yplot);

    % clear NaN
    [fnan1,fnan2] = find(isnan(wx)); wx(fnan1,fnan2)=0;
    [fnan1,fnan2] = find(isnan(wy)); wy(fnan1,fnan2)=0;
    [fnan1,fnan2] = find(isnan(wz)); wz(fnan1,fnan2)=0;

    if multiple
        h1 = quiver3(xplot,yplot,zplot,wx,wy,wz,0);
        lines = findobj(h1,'Type','line');
        for ilines = 1:length(lines)
            set(lines(ilines),'Xdata', get(lines(ilines),'Xdata')/stretch(1));
            set(lines(ilines),'Ydata', get(lines(ilines),'Ydata')/stretch(2));
            set(lines(ilines),'Zdata', get(lines(ilines),'Zdata')/stretch(3));
        end
    else
        h1 = quiver(xplot,yplot, wx,wy,0);
        lines = findobj(h1,'Type','line');
        for ilines = 1:length(lines)
            set(lines(ilines),'Xdata', get(lines(ilines),'Xdata')/stretch(1));
            set(lines(ilines),'Ydata', get(lines(ilines),'Ydata')/stretch(2));
        end
        if max(max(abs(wz))) ~= 0
            h1 = arrow_zplot(xplot,yplot,wz);
            circles = findobj(h1,'Type','line');
            for ipatch = 1:length(circles)
                set(circles(ipatch),'Xdata',get(circles(ipatch),'Xdata')/stretch(1));
                set(circles(ipatch),'Ydata',get(circles(ipatch),'Ydata')/stretch(2));
            end
        end
    end
end
xlim([x1 x2]);
ylim([y1 y2]);
if ratio
    daspect([1 1 1]);
end

if multiple
    izmin = max(1,min(izcut)-1);
    izmax = min(nz,max(izcut)+1);
    zlim([zm(izmin) zm(izmax)]);

    view(40,16)
end

box on
if hhh
    hold on
else
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hh = arrow_zplot(xplot,yplot,z)
    [nxplot,nyplot]=size(xplot);
    xdat = [];
    ydat = [];
    nseg = 20;
    s2 = sqrt(0.5);
    for ix = 1:nxplot
        for iy = 1:nyplot
            rad = abs(z(ix,iy));
            x0 = xplot(ix,iy); y0 = yplot(ix,iy);
            for iseg = 0:nseg
                xdat = [xdat x0+rad*cos(2*pi*iseg/nseg)];
                ydat = [ydat y0+rad*sin(2*pi*iseg/nseg)];
            end
            xdat = [xdat NaN]; ydat = [ydat NaN];
            if z(ix,iy)<0
                xdat = [xdat x0+rad*s2 x0-rad*s2 NaN x0+rad*s2 x0-rad*s2 NaN];
                ydat = [ydat y0+rad*s2 y0-rad*s2 NaN y0-rad*s2 y0+rad*s2 NaN];
            end
        end
    end
    hh = plot(xdat,ydat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
