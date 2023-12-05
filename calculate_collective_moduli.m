function [tvec,Gt,Freq,Storage,Loss,Visc,Comp,tau,X0,Y0] = ...
    calculate_overall_moduli(CondensateChainsG,b,T,phi,xi,t,w)

% Calculates viscoelastic properties of condensate according to
% graph-theoretic formulation of Rouse-Zimm model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of nodes in graph
NumChains = numnodes(CondensateChainsG);

% Graph Laplacian aka Zimm matrix
Zimm = laplacian(CondensateChainsG);

% Find eigenvalues of Zimm matrix
lambda = eig(Zimm);

% Calculate relaxation times, treating simulation T in energy units
tau = xi*b^2*(6*T*lambda(2:end)).^-1; % First eigenvalue is zero

% Calculate viscoelastic properties:

% Zero-shear-rate bulk viscosity in units of pressure x time
Visc = (T*phi/(NumChains*b^3))*sum(tau);

% Steady-state recoverable compliance in units of inverse pressure
Comp = ((NumChains*b^3)/(phi*T))*sum(tau.^2)/sum(tau).^2;

% Calculate the relaxation modulus in units of pressure
Gt(t) = (phi*T/(NumChains*b^3))*sum(exp(-t./tau))*heaviside(t);
% The Heaviside function sets G(t) = 0 for t < 0.

% Take the continuous-time Fourier transform
% (Could probably calculate this analytically to speed up the code.
% See equations in Rouse 1953.)
Gw(w) = fourier(Gt(t),t,w); % Angular frequency, w

% Calculate the average dynamic moduli
Gw1(w) = real(1i*w*Gw(w)); % Storage modulus, G'
Gw2(w) = imag(1i*w*Gw(w)); % Loss modulus, G"

% Convert moduli to type double

% Find range of tau and sample the relaxation modulus
min_order_t = floor(log10(min(tau)));
max_order_t = ceil(log10(max(tau)));
tvec = logspace(min_order_t-1,max_order_t+1,1000); % Covers all tau
Gt = double(Gt(tvec));

% Find range of frequencies and sample the dynamic moduli
min_order_w = floor(log10(2*pi/max(tau)));
max_order_w = ceil(log10(2*pi/min(tau)));
Freq = logspace(min_order_w-3,max_order_w,1000); % All frequencies
Storage = double(Gw1(Freq));
Loss = double(Gw2(Freq));

% Calculate crossover
[X0,Y0] = intersections(Freq,Storage,Freq,Loss,'robust');
