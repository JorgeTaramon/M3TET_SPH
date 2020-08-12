function [iout,xyz_PT] = find_pts_outside_domain(GCOORD,xyz_PT,OPTS)

npt   = size(xyz_PT,2);
iout  = false(1,npt);
xtol  = 1e-12 * (OPTS.xmax-OPTS.xmin);
ytol  = 1e-12 * (OPTS.ymax-OPTS.ymin);

% CHECK IF POINTS ARE INSIDE THE DOMAIN

% xmax
ind = xyz_PT(1,:) > OPTS.xmax;% -xtol;
if any(ind)
    xyz_PT(1,ind) = OPTS.xmax-xtol;
    iout(ind)     = true;
end

% xmin
ind = xyz_PT(1,:) < OPTS.xmin;% +xtol;
if any(ind)
    xyz_PT(1,ind) = OPTS.xmin+xtol;
    iout(ind)     = true;
end

% ymax
ind = xyz_PT(2,:) > OPTS.ymax;% -ytol;
if any(ind)
    xyz_PT(2,ind) = OPTS.ymax-ytol;
    iout(ind)     = true;
end

% ymin
ind = xyz_PT(2,:) < OPTS.ymin;% +ytol;
if any(ind)
    xyz_PT(2,ind) = OPTS.ymin+ytol;
    iout(ind)     = true;
end

% zmax
if isfield(OPTS,'EL2NOD_top') && ~isempty(OPTS.EL2NOD_top)
    zmax = zcoord_boundary(GCOORD,OPTS.EL2NOD_top,xyz_PT);
    ind  = xyz_PT(3,:) > zmax-OPTS.ztol;
    if any(ind)
        % Adjust all that are very close to zmax but...
        xyz_PT(3,ind) = zmax(ind);
        % ...only mark the ones that are really above top boundary
        ind           = xyz_PT(3,:) > zmax;
        iout(ind)     = true;
    end
else
    ind  = xyz_PT(3,:) > OPTS.zmax;
    if any(ind)
        xyz_PT(3,ind) = OPTS.zmax;
        iout(ind)     = true;
    end
end

% zmin
if  isfield(OPTS,'EL2NOD_bot') && ~isempty(OPTS.EL2NOD_bot)
    zmin = zcoord_boundary(GCOORD,OPTS.EL2NOD_bot,xyz_PT);
    ind  = xyz_PT(3,:) < zmin+OPTS.ztol;
    if any(ind)
        % Adjust all that are very close to zmin but...
        xyz_PT(3,ind) = zmin(ind);
        % ...only mark the ones that are really below bottom boundary
        ind           = xyz_PT(3,:) < zmin; 
        iout(ind)     = true;
    end
else
    ind  = xyz_PT(3,:) < OPTS.zmin;
    if any(ind)
        xyz_PT(3,ind) = OPTS.zmin;
        iout(ind)     = true;
    end
end

end % END OF FUNCTION find_pts_outside_domain