function NUMSCALE = numerical_eq_scaling(scaling)

%==========================================================================
% SCALING OF EQUATIONS TO REDUCE NUMERICAL ROUND-OFF PROBLEMS
%
% This is needed if not all input parameters are in SI-units (e.g. length
% in km, time in Myr, viscosity in multiples of 1e19 Pa-s,...).
% Scaling of the viscous flow equations results in a single term "Bscale"
% that multiplies the buoyancy terms in the equations. Thermal diffusivity
% is also scaled so that it matches the time and length units t0 and L0, 
% respectively. 
%
% JH Jul 2014

if nargin==0
    error('Define scaling method as input argument.');
end

switch scaling
    case 'mantle'
        Myr             = 1e6*365.25*24*60*60;
        NUMSCALE.L0     = 1000;  % m; length unit
        NUMSCALE.t0     = Myr;   % s; time unit
        NUMSCALE.Visc0  = 2e18;  %1e19; %1e17;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 3300;  % 3428.7;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 'Myr';
        NUMSCALE.unit_L = 'km';
        
    case 'earth'
        Myr             = 1e6*365.25*24*60*60;
        NUMSCALE.L0     = 1e6;   % m; length unit
        NUMSCALE.t0     = Myr;   % s; time unit
        NUMSCALE.Visc0  = 1e21;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 3300;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 'Myr';
        NUMSCALE.unit_L = 'x 1000 km';
        
    case 'SI' % SI units
        NUMSCALE.L0     = 1;  % m; length unit
        NUMSCALE.t0     = 1;  % s; time unit
        NUMSCALE.Visc0  = 1;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 0;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 's';
        NUMSCALE.unit_L = 'm';

    otherwise
        error(' Unknown scaling method. See function "numerical_eq_scaling" for available methods.');
end

%==========================================================================
% Do not edit the next lines...

NUMSCALE.U0     = NUMSCALE.L0 / NUMSCALE.t0;
    % m/s; velocity unit
NUMSCALE.unit_U = [NUMSCALE.unit_L '/' NUMSCALE.unit_t];

NUMSCALE.Kappa0 = NUMSCALE.L0^2 / NUMSCALE.t0;
    % calculate diffusivity unit
NUMSCALE.P0     = NUMSCALE.Visc0 / NUMSCALE.t0;
    % scaling for pressure and shear modulus

% Calculate buoyancy scaling factor
% This factor multiplies the gravitaional force in the element assembly
NUMSCALE.Bscale = NUMSCALE.L0^2 / (NUMSCALE.Visc0 * NUMSCALE.U0);

%==========================================================================

end % END OF FUNCTION numerical_eq_scaling