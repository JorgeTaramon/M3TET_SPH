function [els,lc] = locate_points_in_neighboring_element(els_wrng,gTH_PT,GCOORD_SPH,EL2NOD,MESH)

tol_step  = 1e-6;
npnt_wrng = length(els_wrng);
els       = zeros(1,npnt_wrng);
lc        = zeros(3,npnt_wrng);
fprintf(' %1i points have not been located by tsearch2\n',npnt_wrng);
for ipnt=1:npnt_wrng
    iel     = els_wrng(ipnt);
    fprintf(' Point %3i...',ipnt);
    if iel>0
        % Find all elements connected to the wrong element
        el_nb   = find(any(ismember(EL2NOD,EL2NOD(1:4,iel)),1));
        nel_nb  = length(el_nb);

        % Assume the point is in all the element neighbors and calculate local coordinates
        lc_nb      = ...
            local_coords_3d_sph(GCOORD_SPH,EL2NOD,el_nb,repmat(gTH_PT(:,ipnt),1,nel_nb),MESH);
        lc_nb(4,:) = 1-sum(lc_nb(1:3,:));
        iel_ok_nb  = find(all(lc_nb >= 0 & lc_nb <= 1,1));
        if isempty(iel_ok_nb)
            % try to find the right element using a tolerance for the lc 
            tol = 0;
            while isempty(iel_ok_nb) && tol <= 15e-3
                tol       = tol + tol_step;
                iel_ok_nb = find(all(lc_nb >= 0-tol & lc_nb <= 1+tol,1));
            end
        end
    else
        iel_ok_nb = [];
    end
    % If the point was not located in the neighboring elements we have to search in all elements
    if isempty(iel_ok_nb)
        lc_all        = ...
            local_coords_3d_sph(GCOORD_SPH,EL2NOD,1:MESH.nel,repmat(gTH_PT(:,ipnt),1,MESH.nel),MESH);
        lc_all(4,:)   = 1-sum(lc_all(1:3,:));
        iel_ok_non_nb = find(all(lc_all >= 0 & lc_all <= 1,1));
        if isempty(iel_ok_non_nb)
            tol = 0;
            % try to find the right element using a tolerance for the lc 
            while isempty(iel_ok_non_nb) && tol <= 15e-3
                tol           = tol + tol_step;
                iel_ok_non_nb = find(all(lc_all >= 0-tol & lc_all <= 1+tol,1));
            end
        end
        if isempty(iel_ok_non_nb)
            error(' This should not happen');
        elseif length(iel_ok_non_nb)>1
            if sum(ismember(iel_ok_non_nb,MESH.els_in_cone_iso)) == length(iel_ok_non_nb) % -> all elements are isoparametric
                % In order to decide which iel_ok is the right one, we need to check which one has the smallest delS
                gTH_PT_els_in_cone_iso     = gTH_PT(:,ipnt); % coordinates of points to be located
                gX_PT_els_in_cone_iso      = spherical2cartesian(gTH_PT_els_in_cone_iso);
                RR_X_90_CCW                = [ 1  0  0 ; 0  0 -1 ; 0  1  0 ];                % rotation matrix
                gX_PT_els_in_cone_iso_rot  = RR_X_90_CCW * gX_PT_els_in_cone_iso;            % rotate those points
                gTH_PT_els_in_cone_iso_rot = cartesian2spherical(gX_PT_els_in_cone_iso_rot); % convert back to spherical coordinates
                GCOORD_rot_90_X            = RR_X_90_CCW * MESH.GCOORD;                      % rotate all mesh nodes in Cartesian coordinates
                GCOORD_SPH_rot_90_X        = cartesian2spherical(GCOORD_rot_90_X);           % convert back to spherical coordinates
                [~,delS]                   = ...
                    local_coords_curved_3d_sph(GCOORD_SPH_rot_90_X,EL2NOD,iel_ok_non_nb,repmat(gTH_PT_els_in_cone_iso_rot,1,length(iel_ok_non_nb)));
                iel_ok_non_nb              = iel_ok_non_nb(sum(abs(delS))==min(sum(abs(delS))));
                fprintf('located in non-neighboring element\n');
                % Save correct element
                els(ipnt)  = iel_ok_non_nb;
                lc(:,ipnt) = lc_all(1:3,iel_ok_non_nb);
            else % 2 or more elements (isoparametric and non-isoparametric)
                iel_ok_non_nb(ismember(iel_ok_non_nb,MESH.els_in_cone_iso)) = []; % Remove from the list isoparametric elements since they have curved edges (less precision)
                if length(iel_ok_non_nb) == 1
                    fprintf('located in non-neighboring element\n');
                    % Save correct element
                    els(ipnt)  = iel_ok_non_nb;
                    lc(:,ipnt) = lc_all(1:3,iel_ok_non_nb);
                elseif length(iel_ok_non_nb) == 2 % check if the elements share edges, if so it means that the point is located on one of the shared edges
                    if sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(1))),GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(2))))) == 3) > 0
                        fprintf('located in non-neighboring element\n');
                        % Save correct element
                        els(ipnt)  = iel_ok_non_nb(1); % pick one of them
                        lc(:,ipnt) = lc_all(1:3,iel_ok_non_nb(1));
                    else
                        error('located in %1i separated elements',length(iel_ok_non_nb));
                    end
                elseif length(iel_ok_non_nb) == 3 % check if the elements share edges, if so it means that the point is located on one of the shared edges
                    if sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(1))),GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(2))))) == 3) > 0 && ...
                       sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(2))),GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(3))))) == 3) > 0 && ...
                       sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(3))),GCOORD_SPH(:,EL2NOD(1:4,iel_ok_non_nb(1))))) == 3) > 0
                        fprintf('located in non-neighboring element\n');
                        % Save correct element
                        els(ipnt)  = iel_ok_non_nb(1); % pick one of them
                        lc(:,ipnt) = lc_all(1:3,iel_ok_non_nb(1));
                    else
                        error('located in %1i separated elements',length(iel_ok_nb));
                    end
                else
                    error('located in %1i elements',length(iel_ok_non_nb));
                end
            end
        else
            fprintf('located in non-neighboring element\n');
            % Save correct element
            els(ipnt)  = iel_ok_non_nb;
            lc(:,ipnt) = lc_all(1:3,iel_ok_non_nb);
        end
    elseif length(iel_ok_nb)>1
        if sum(ismember(el_nb(iel_ok_nb),MESH.els_in_cone_iso)) == length(iel_ok_nb) % -> all elements are isoparametric
            % In order to decide which iel_ok is the right one, we need to check which one has the smallest delS
            gTH_PT_els_in_cone_iso     = gTH_PT(:,ipnt); % coordinates of points to be located
            gX_PT_els_in_cone_iso      = spherical2cartesian(gTH_PT_els_in_cone_iso);
            RR_X_90_CCW                = [ 1  0  0 ; 0  0 -1 ; 0  1  0 ];                % rotation matrix
            gX_PT_els_in_cone_iso_rot  = RR_X_90_CCW * gX_PT_els_in_cone_iso;            % rotate those points
            gTH_PT_els_in_cone_iso_rot = cartesian2spherical(gX_PT_els_in_cone_iso_rot); % convert back to spherical coordinates
            GCOORD_rot_90_X            = RR_X_90_CCW * MESH.GCOORD;                      % rotate all mesh nodes in Cartesian coordinates
            GCOORD_SPH_rot_90_X        = cartesian2spherical(GCOORD_rot_90_X);           % convert back to spherical coordinates
            [~,delS]                   = ...
                local_coords_curved_3d_sph(GCOORD_SPH_rot_90_X,EL2NOD,el_nb(iel_ok_nb),repmat(gTH_PT_els_in_cone_iso_rot,1,length(iel_ok_nb)));
            iel_ok_nb                  = iel_ok_nb(sum(abs(delS))==min(sum(abs(delS))));
            fprintf('located in neighboring element\n');
            % Save correct element
            els(ipnt)  = el_nb(iel_ok_nb);
            lc(:,ipnt) = lc_nb(1:3,iel_ok_nb);
        else % 2 or more elements (isoparametric and non-isoparametric)
            iel_ok_nb(ismember(iel_ok_nb,MESH.els_in_cone_iso)) = []; % Remove from the list isoparametric elements since they have curved edges (less precision)
            if length(iel_ok_nb) == 1
                fprintf('located in neighboring element\n');
                % Save correct element
                els(ipnt)  = el_nb(iel_ok_nb);
                lc(:,ipnt) = lc_nb(1:3,iel_ok_nb);
            elseif length(iel_ok_nb) == 2 % check if the elements share edges, if so it means that the point is located on one of the shared edges
                if sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(1)))),GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(2)))))) == 3) > 0
                    fprintf('located in neighboring element\n');
                    % Save correct element
                    els(ipnt)  = el_nb(iel_ok_nb(1)); % pick one of them
                    lc(:,ipnt) = lc_nb(1:3,iel_ok_nb(1));
                else
                    error('located in %1i separated elements',length(iel_ok_nb));
                end
            elseif length(iel_ok_nb) == 3 % check if the elements share edges, if so it means that the point is located on one of the shared edges
                if sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(1)))),GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(2)))))) == 3) > 0 && ...
                   sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(2)))),GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(3)))))) == 3) > 0 && ...
                   sum(sum(ismember(GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(3)))),GCOORD_SPH(:,EL2NOD(1:4,el_nb(iel_ok_nb(1)))))) == 3) > 0
                    fprintf('located in neighboring element\n');
                    % Save correct element
                    els(ipnt)  = el_nb(iel_ok_nb(1)); % pick one of them
                    lc(:,ipnt) = lc_nb(1:3,iel_ok_nb(1));
                else
                    error('located in %1i separated elements',length(iel_ok_nb));
                end
            else
                error('located in %1i elements',length(iel_ok_nb));
            end
        end
    else
        fprintf('located in neighboring element\n');
        % Save correct element
        els(ipnt)  = el_nb(iel_ok_nb);
        lc(:,ipnt) = lc_nb(1:3,iel_ok_nb);
    end    
end

end