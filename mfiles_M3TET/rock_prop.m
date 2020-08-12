function VAR = rock_prop(VAR,PHYSICS,NUMSCALE,MESH,SETTINGS)
% Usage: VAR = rock_prop(VAR,PHYSICS,NUMSCALE,MESH,SETTINGS)
% 
% Purpose: Update rock density and viscosity.
%
% Input:
%   VAR      : [structure] : major variable fields, each is a vector
%   PHYSICS  : [structure] : physical properties
%   SETTINGS : [structure] : model parameters
%   NUMSCALE : [structure] : numerical scaling parameters
%   MESH     : [structure] : FE mesh parameters
%
% Output:
%   VAR      : [structure] : major variable fields, each is a vector
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JH Jul 2013: introduced T_K and T0_K (VAR.T and T0 in units of Kelvin)
% JH Aug 2014: added element Phases (i.e. different rock types)
%

fid_log = SETTINGS.fid_log;
fprintf(fid_log,'\n Updating rock properties...');

% Reference values
% =================
Visc0  = PHYSICS.Visc;
Dens0  = PHYSICS.Dens;
T0     = PHYSICS.T_ref;
radius = sqrt( sum(MESH.GCOORD.^2,1) );
depth  = max(radius)-radius;

% =========================================================================
% DENSITY
% =========================================================================

% Temperature effect on density
% =============================
% (alpha is thermal expansion coeff, VAR.Dens is defined at the nodes)
VAR.Dens = Dens0 * (1 - PHYSICS.alpha*(VAR.T-T0)); % new density

% Density linearly varying with depth for the top half and constant for the bottom half 
% a = (Dens0 - Dens0*(1.1))/(6371-4921);
% b = Dens0 - a*6371;
% VAR.Dens = a*radius' + b;
% VAR.Dens(radius<4921) = 1.1*Dens0;

% Density cut-offs (if defined)
% =============================
if isfield(PHYSICS,'maxDens') && ~isempty(PHYSICS.maxDens)
    VAR.Dens = min(VAR.Dens,PHYSICS.maxDens); % upper density cut-off
end
if isfield(PHYSICS,'minDens') && ~isempty(PHYSICS.minDens)
    VAR.Dens = max(VAR.Dens,PHYSICS.minDens); % lower density cut-off
end

% =========================================================================
% VISCOSITY
% =========================================================================
vALL = ones(MESH.nnod,1);
R    = 8.314472; % J/mol K; univ. gas constant

% Temperature effect on viscosity
% ===============================
T0_K = T0 + 273;   % T0 in Kelvin
Ea   = PHYSICS.Ea; % J/mol; activation energy
if Ea>0
    if Ea<1000
        error(' Activation energy PHYSICS.Ea seems to be in kJ/mol. Convert to J/mol.');
    end
    T_K  = VAR.T + 273;   % T in Kelvin
%     T_K  = T_K + 0.3*z;   % adiabate
    eta  = Ea/(R*T0_K);
    vALL = vALL .* exp( eta*( (T0_K./T_K) - 1) );
        % T-effect on viscosity; see Turcotte & Schubert, Geodynamics
        % Factor by which viscosity changes at every node
    vALL = min(vALL,1e8);
    vALL = max(vALL,1e-8);
        % limit viscosity change
end

% Pressure (depth) effect on viscosity (see Hirth & Kohlstedt 2003)
% =================================================================
vALL2 = vALL;
Va    = PHYSICS.Va; % m3/mol; activation volume
if Va>0
    P = abs(PHYSICS.g * Dens0 * NUMSCALE.L0 * depth(:));
        % lithostatic pressure at z
    vALL2 = vALL2 .* exp( P.*Va ./ (R.*T0_K) );
        % P-effect on viscosity (factor)
end

% Viscosity increase in lower mantle
% ==================================
vALL3 = vALL;
if isfield(PHYSICS,'a_UM2LM') && PHYSICS.a_UM2LM>0
    vALL3 = vALL3 * PHYSICS.a_UM2LM;
end

% Compute the minimum between pressure (depth) effect on viscosity and viscosity increase in the lower mantle 
% ===========================================================================================================
vALL = min([vALL2 vALL3],[],2); % This is done to avoid the viscosity "jump" in the upper mantle - lowe mantle transition

% New NODAL viscosity
% ===================
VAR.Visc = vALL .* Visc0;


% Viscosity cut-offs (if defined)
% ===============================
if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
    VAR.Visc = min(VAR.Visc,PHYSICS.maxVisc); % upper viscosity cut-off
end
if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
    VAR.Visc = max(VAR.Visc,PHYSICS.minVisc); % lower viscosity cut-off
end

fprintf(fid_log,' done.\n\n');

end % END OF FUNCTION rock_prop