function [tvec,Gt,Freq,Storage,Loss,avg_Visc,std_Visc,avg_Comp,...
    std_Comp,tau_all,std_tau_all,X0,Y0,ChainsSubG] = ...
    calculate_single_chain_moduli(CondensateSitesG,ind_chain,Sites,...
    siteType,b,xi,T,phi,RelEl,t,w)

% Calculates viscoelastic moduli for each chain according to
% graph-theoretic formulation of Rouse-Zimm model. Uses edge-weighted
% single-chain sub-graphs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumChains = length(ind_chain); % Number of chains in condensate

tau = cell(NumChains,1);
Visc = zeros(NumChains,1);
Comp = zeros(NumChains,1);

for ii = 1:NumChains

    % Find sites associated with each chain in condensate:
    AllSitesInChain = linspace(Sites*ii-Sites+1,Sites*ii,Sites);
    % Generate sub-graph from topology of full graph:
    ChainsSubG = subgraph(CondensateSitesG,AllSitesInChain);
    % Edges in sub-graph represent intrachain interactions only.
    
    % Weight edges by absolute value of relative energies:

    % Find edges
    edgeID = table2array(ChainsSubG.Edges);

    % Find corresponding site types
    edge_type = siteType(edgeID);
    edge_type(edge_type == 0) = 11; % For A1-LCD only, which has 1-11 types of residues
    % Sort each row in ascending order
    edge_type = sort(edge_type,2);

    % Assign edge weights
    idx = sub2ind(size(RelEl),edge_type(:,1),edge_type(:,2)); % Index
    edge_energies = RelEl(idx); % Energy associated with each edge
    W = diag(abs(edge_energies)); % Matrix of edge weights
    
    % Calculate graph laplacian, aka Zimm matrix, for weighted graph
    B = incidence(ChainsSubG); % Unweighted incidence matrix
    Zimm = B*W*B'; % Zimm matrix

    % Find eigenvalues of Zimm matrix
    lambda = eig(Zimm);

    % Calculate relaxation times, treating simulation T in energy units
    dummy_tau = xi*b^2*(6*T*lambda(2:end)).^-1; % First eigenvalue is zero
    tau{ii} = dummy_tau;

    % Calculate viscoleastic properties:
    
    % Zero-shear-rate bulk viscosity in units of pressure x time
    Visc(ii) = (T*phi/(NumChains*b^3))*sum(dummy_tau);
    
    % Steady-state recoverable compliance in units of inverse pressure
    Comp(ii) = ((NumChains*b^3)/(phi*T))*sum(dummy_tau.^2)/sum(dummy_tau).^2;

end

% Average time constants over all chains
tau_all = mean([tau{:}],2);
% Standard deviation
std_tau_all = std([tau{:}],0,2);

% Calculate the relaxation modulus in units of pressure
Gt(t) = (phi*T/(NumChains*b^3))*sum(exp(-t./tau_all))*heaviside(t);
% The Heaviside function sets G(t) = 0 for t < 0.

% Take the continuous-time Fourier transform
Gw(w) = fourier(Gt(t),t,w); % Angular frequency w

% Calculate the average dynamic moduli
Gw1(w) = real(1i*w*Gw(w)); % Storage modulus, G'
Gw2(w) = imag(1i*w*Gw(w)); % Loss modulus, G"

% Calculate average viscosity per snapshot
avg_Visc = mean(Visc);
% Calculate standard deviation per snapshot
std_Visc = std(Visc);

% Calculate average compliance per snapshot
avg_Comp = mean(Comp);
% Calculate standard deviation per snapshot
std_Comp = std(Comp);

% Convert moduli to type double

% Find range of tau and sample the relaxation modulus
min_order_t = floor(log10(min(tau_all)));
max_order_t = ceil(log10(max(tau_all)));
tvec = logspace(min_order_t-1,max_order_t+1,1000); % Covers all tau
Gt = double(Gt(tvec));

% Find range of frequencies and sample the dynamic moduli
min_order_w = floor(log10(2*pi/max(tau_all)));
max_order_w = ceil(log10(2*pi/min(tau_all)));
Freq = logspace(min_order_w-3,max_order_w,1000); % All frequencies
Storage = double(Gw1(Freq));
Loss = double(Gw2(Freq));

% Calculate crossover
[X0,Y0] = intersections(Freq,Storage,Freq,Loss,'robust');
