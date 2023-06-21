function [CondensateSitesG,ind_chain,siteType] = ...
    generate_site_graph(A,Length,Sites)

% Makes graph of all sites in condensate with respect to non-bonded
% interactions, plus chain connectivity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find which chains are in the condensate by first making graph of all
% chains in the system and then determining largest connected network. By
% system, we mean dilute + dense phase.

% Calculate distances between each site

SitesInSystem = length(A); % Number of sites in system

Ai = reshape(A(:,4:6),[],1,3);
B = zeros(SitesInSystem,SitesInSystem);
for j = 1:SitesInSystem
    Aj = reshape(A(j,4:6),1,[],3);
    dxyz = Ai-Aj; % Differences dx, dy, dz
    % Apply minimum image convention
    dxyz = mod(dxyz+Length/2,Length)-Length/2;
    B(:,j) = sqrt(sum(dxyz.^2,3));
end

% Generate graph of all chains in the system: nodes = chains

% Calculate adjacency matrix for sites. Do this by applying distance
% criteria to find adjacent sites and then converting to 1's and 0's:
C = ((B < 1.75) & (B > 0));
% Adjacency matrix for chains: N x N matrix s.t. nij = 1 if two chains are
% adjacent if any sites between them are adjacent and 0 otherwise:
L = A(:,3)+1; % Chain ID for each site
S = sparse(L,1:numel(L),1); % sij = 1 if site is in chain, 0 otherwise
D = S*C*S' > 0; % Adjacency matrix except dii = 1
D(1:size(D,1)+1:end) = 0; % Adjacency matrix for chains
SystemChainsG = graph(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate graph of all chains in condensate: nodes = chains

% An edge in SystemChainsG means that at least one site on two chains is
% adjacent. Use something like union-find to find the largest sub-graph.
bincell = biconncomp(SystemChainsG, 'OutputForm', 'cell');
[s,d] = cellfun(@size,bincell);
out = max([s,d],[],'omitnan');
ind = d == out | s == out;
bincell_largest = bincell(ind);
ind_chain = cell2mat(bincell_largest); % Which nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate adjacency matrix for single chain in condensate: nodes = sites

% Determine all sites in condensate
SitesInCondensate = find(ismember(A(:,3)+1,ind_chain));

% Modify graph such that the edges represent covalent bonding within each
% chain as well as non-bonded interactions. Treat all of the interactions
% depicted by the graph as being the same.

% Build adjacency matrix with edges as covalent bonds between sites. Note
% that in the LaSSI set up, covalent bonding is specified in sequential
% order for each site ID for each chain in the trajectory file.

% Assign edges
s = linspace(1,length(SitesInCondensate)-1,length(SitesInCondensate)-1);
s(Sites:Sites:end) = [];

t = linspace(2,length(SitesInCondensate),length(SitesInCondensate)-1);
t(Sites:Sites:end) = [];

% Build graph and take adjacency matrix
P = graph(s,t);
C_bond = adjacency(P);

% Separately, build adjacency matrix for non-bonded interactions in
% condensate. Here, number of nodes = number of sites in the dense phase.

% Get adjacency matrix for sites in condensate
SitesInDilute = find(~ismember(A(:,1)+1,SitesInCondensate)); % with respect to sites IDs
% Remove dilute-phase sites from the site-site adjacency matrix
C(ismember(1:size(C,1),SitesInDilute),:) = [];
C(:,ismember(1:size(C,2),SitesInDilute)) = [];
% Order sites wrt condensate
siteIDs = linspace(1,length(SitesInCondensate),length(SitesInCondensate));
C_nonbonded = sparse(length(SitesInCondensate),length(SitesInCondensate));
C_nonbonded(siteIDs,siteIDs) = C_nonbonded(siteIDs,siteIDs) | C;

% Combine adjacency matrices
C_condensate = (C_bond | C_nonbonded == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate graph for all sites in condensate, where edges represent
% non-bonded interactions and chain connectivity
CondensateSitesG = graph(C_condensate);

% Remove self-loops
CondensateSitesG = simplify(CondensateSitesG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find site type for each site on single chain

siteType = A(1:Sites,2); % Will need to weight the sub-graphs
