function [p, uy] = terzaghiReference(z, t, H, q, E, nu, k, mu, alpha, Sm, np)
% Analytical solution for a homogeneous, laterally confined poroelastic
% layer.
% Inputs:
%   z: elevation above the fixed, impermeable base (m), 0 <= z <= H.
%   t: elapsed time since an instantaneous, constant compressive load (s).
%   H: layer thickness (m); 
%   q: applied compressive traction magnitude (Pa).
%   E: Young's modulus (Pa);
%   nu: Poisson's ratio;
%   k: permeability (m2).
%   mu: fluid viscosity (Pa*s);
%   alpha: Biot coefficient;
%   Sm: storage (1/Pa).
%   np: number of Fourier terms (default 1000).
% Outputs:
%   Column vectors: excess pressure p (Pa) and displacement uy (m).
% Conventions:
%   - Displacement is positive upwards.
%   - The top is drained.
%   - Gravity is neglected.
%   - Solid skeleton is linear elastic with zero lateral strains.
%
% Reference: Verruijt, A. (2016), Theory and Problems of Poroelasticity,
% Section 2.2, Eqs. (2.2), (2.6), (2.9), (2.10), and (2.25).
% Includes compressible fluid and grains through alpha and Sm.
% Here H = h, Sm = S, and k/mu corresponds to the book's k/gamma_f.
% As in the book, z is measured upwards from the base.
% Displacement follows from Eq. (2.37):
% uy_z = mv*(alpha*p-q), with uy(0,t) = 0.
% The book's positive surface settlement is w = -uy(H,t).

if nargin < 11
    np = 1000;
end

validateattributes(H, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(z, {'numeric'}, {'vector','real','finite','nonempty','>=',0,'<=',H});
validateattributes(t, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(q, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(E, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(nu, {'numeric'}, {'scalar','real','finite','>',-1,'<',0.5});
validateattributes(k, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(mu, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(alpha, {'numeric'}, {'scalar','real','finite','nonnegative','<=',1});
validateattributes(Sm, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(np, {'numeric'}, {'scalar','real','finite','integer','positive'});

% Confined compressibility, consolidation coefficient and initial pressure
mv = (1+nu)*(1-2*nu)/(E*(1-nu));
storage = Sm + alpha*alpha*mv;
assert(storage > 0, 'The combined storage coefficient must be positive.');
Cv = k/(mu*storage);
p0 = alpha*mv*q/storage;
z = z(:);

% The initial pressure is uniform inside the layer, with a drained top.
if t == 0
    p = p0*ones(size(z));
    p(z == H) = 0;
    uy = mv*(alpha*p0-q)*z;
    return;
end

% Pressure series and its integral from the base to each elevation
p = zeros(size(z));
integralPressure = zeros(size(z));
T = Cv*t/(H*H);
for j = 1:np
    M = (2*j - 1.0) * pi/2.0;
    coefficient = 2.0 * p0 * (-1)^(j-1)/M*exp(-M*M*T);
    p = p + coefficient*cos(M*z/H);
    integralPressure = integralPressure + coefficient*H/M*sin(M*z/H);
end
p(z == H) = 0;
uy = mv*(alpha*integralPressure-q*z);
end
