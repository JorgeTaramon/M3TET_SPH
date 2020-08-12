function V_i = interp3d_tet(EL2NOD,els,lc,V,OPTS)
% Usage: V_i = interp3d_tet(EL2NOD,els,lc,V,OPTS)
%
% Purpose: Interpolates variables in a triangular finite element mesh using
%          linear or quadratic interpolation functions.
%
% Input
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   els    : [rowvector] : elements in which will be interpolated
%   lc     : [matrix]    : local coordinates in each element
%   V      : [matrix]    : variable field(s) to be interpolated (nnod x nvar)
%   OPTS   : options
%
% Output
%   V_i    : [colvector] : interpolated values
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JH Dec 2012 : function will interpolate linearly (4-node tetrahedra) and
%               quadratically (10-node tetrahedra) depending on size
%               of EL2NOD
% JH Jan 2013 : handles NaN in lc
% JH Jan 2014 : uses shape functions from sf_dsf_tet410
%               (4- and 10-node intepolation possible)
%
                                                                            t0 = tic;
% DEFAULT SETTINGS (will be overwritten if options (OPTS) are provided
if nargin==4; OPTS = []; end
if ~isfield(OPTS,'verbose'); OPTS.verbose = 0; end
if ~isfield(OPTS,'monotonic'); OPTS.monotonic = 0; end


ind  = ~isnan(lc(1,:)) & ~isnan(lc(2,:)) & ~isnan(lc(3,:));
els  = els(ind);
if iscell(EL2NOD); EL2NOD = EL2NOD{1}; end
nnod = max(max(EL2NOD));
nnodel = size(EL2NOD,1);

transpose_result = 0;
if size(V,2)==nnod
    V=V'; transpose_result = 1;
end
nvar = size(V,2);
npt  = length(els);

% Local coordinates
N    = sf_dsf_tet(lc(:,ind),nnodel,'matrix');

% Interpolate to points using shape functions
V_i  = nan(npt,nvar);
for ivar=1:nvar
    B   = V(:,ivar);
    B_i = sum(N.*B(EL2NOD(:,els)))';
    if nnodel==10 && OPTS.monotonic
        B_elnods = B(EL2NOD(1:nnodel,:));
        B_max_el = max(B_elnods,[],1)';
        B_min_el = min(B_elnods,[],1)';
        B_max    = B_max_el(els);
        B_min    = B_min_el(els);
        B_i      = max(B_i,B_min);
        B_i      = min(B_i,B_max);
    end
    V_i(:,ivar) = B_i(:);
end

if transpose_result
    V_i = V_i';
end

if OPTS.verbose
    fprintf(' Interpolation using %1i shape functions (monotonic=%1i):\n',...
            nnodel,OPTS.monotonic);
    fprintf('      %1i variable(s) at %1i points in %5.2e sec\n',...
                 nvar,npt,toc(t0));
end

end % END OF FUNCTION interp3d_tet