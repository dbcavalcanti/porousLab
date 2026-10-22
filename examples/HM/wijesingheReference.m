function [b_mm, p] = wijesingheReference(x, t, mu, Kn, bi, pi, p0)
% Semi-analytical solution for a semi-infinite fracture with no matrix leak-off.
% x: nonnegative distances from the inlet (m);
% t: time since the pressure step (s), t > 0.
% mu: viscosity (Pa*s);
% Kn: normal stiffness (Pa/m);
% bi: initial aperture (m).
% pi and p0: initial and inlet pressures (Pa).
% Outputs are column vectors: aperture b_mm (mm) and pressure p (Pa).
%
% Reference: A. M. Wijesinghe, "A Similarity Solution for Coupled Deformation
% and Fluid Flow in Discrete Fractures", Eqs. (6)-(17), pp. 351-352.
%
% b_t = Kn/(12*mu)*(b^3*b_x)_x;
% p-pi = Kn*(b-bi).
% eta=x/sqrt(D*t),
% D=bi^3*Kn/(24*mu),
% s=(b/bi)^4:
% s'' = -(eta/4)*s^(-3/4)*s',
% s(0)=(b(0)/bi)^4, s(infinity)=1.

% The similarity coordinate is defined only for positive time.
validateattributes(t, {'numeric'}, {'scalar','real','finite','positive'});

% Inlet aperture ratio b0/bi, obtained from the normal stiffness relation.
w0 = 1+(p0-pi)/(Kn*bi);
assert(w0 > 0, 'The inlet aperture must be positive.');

% Similarity scaling used in the paper: eta = x/sqrt(D*t).
D = bi^3*Kn/(24*mu);

% Approximate the semi-infinite domain by a finite interval in eta.
% The initial mesh can be refined by bvp4c during the solution.
etaMax = 24; % checked against 36 for this benchmark (w0=1.9)
etaMesh = linspace(0,etaMax,1200);
guess = @(eta) wijesingheGuess(eta,w0);

% Solve for y = [s; ds/deta], where s = (b/bi)^4.
% The eps guard avoids nonpositive trial values in the fractional power.
ode = @(eta,y) [y(2); -eta/4*max(y(1),eps)^(-3/4)*y(2)];

% Prescribe the inlet aperture and the initial aperture at the far end.
bc = @(left,right) [left(1)-w0^4; right(1)-1];
opts = bvpset('RelTol',1e-8,'AbsTol',1e-11,'NMax',10000);
sol = bvp4c(ode,bc,bvpinit(etaMesh,guess),opts);

% Evaluate the similarity profile at the requested distances and time.
% Beyond etaMax, use the undisturbed state s = 1 (b = bi and p = pi).
eta = x(:)/sqrt(D*t);
s = ones(size(eta));
inside = eta <= etaMax;
if any(inside)
    s(inside) = deval(sol,eta(inside),1);
end
assert(all(isfinite(s) & s > 0), 'Nonphysical reference aperture.');

% Recover aperture and pressure from s. Keep b in metres for the stiffness
% relation; only the aperture output is converted to millimetres.
b = bi*s.^0.25;
b_mm = b*1e3;
p = pi+Kn*(b-bi);
end

function y = wijesingheGuess(eta,w0)
% Initial guess for the nonlinear boundary value problem.

% Use the small-deformation profile and its derivative with respect to eta.
% This erfc profile is only used to initialize the nonlinear solve.
w = 1+(w0-1)*erfc(eta/sqrt(8));
dw = (1-w0)*exp(-eta.^2/8)/sqrt(2*pi);

% Transform to the solver variables: s = w^4 and ds/deta = 4*w^3*dw/deta.
y = [w.^4; 4*w.^3.*dw];
end
