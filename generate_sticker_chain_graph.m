function CondensateChainsG = generate_sticker_chain_graph(A,Length,Stickers)

% Makes graph of all chains in condensate with respect to stickers only

% Takes as input A = data in LAMMPS trajectory format, Length = length of
% cubic box on each side, Stickers = vector of sticker IDs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find which chains are in the condensate by first making graph of all
% chains in the system and then determining largest connected network. By
% system, we mean dilute + dense phase.

% Calculate distances between each site

% Analyze stickers only. Treat like and unlike interactions as the same.
A = A(ismember(A(:,2),Stickers),:);

SitesInSystem = length(A); % Number of stickers in system

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

% An edge in SystemG means that at least one site on two chains is
% adjacent. Use something like union-find to find the largest sub-graph.
bincell = biconncomp(SystemChainsG, 'OutputForm', 'cell');
[s,d] = cellfun(@size,bincell);
out = max([s,d],[],'omitnan');
ind = d == out | s == out;
bincell_largest = bincell(ind);
% Generate largest sub-graph. Will be unique here.
CondensateChainsG = subgraph(SystemChainsG,cell2mat(bincell_largest));