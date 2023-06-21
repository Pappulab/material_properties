function [CondensateSitesG,ind_chain,siteType] = ...
    generate_sticker_site_graph(A,Length,Sites,Stickers)

% Makes graph of all sites in condensate with respect to non-bonded
% sticker-sticker interactions, plus chain connectivity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find which chains are in the condensate by first making graph of all
% chains in the system and then determining largest connected network. By
% system, we mean dilute + dense phase.

% Calculate distances between each site

% Will need all sites later
A_all = A;

% Analyze stickers only. Treat like and unlike interactions as the same.
StickSitesInSystem = A(ismember(A(:,2),Stickers),:);

N = length(StickSitesInSystem); % Number of stickers in system
Ai = reshape(StickSitesInSystem(:,4:6),[],1,3);
B = zeros(N,N);

for j = 1:N
    Aj = reshape(StickSitesInSystem(j,4:6),1,[],3);
    dxyz = Ai-Aj; % Differences dx, dy, dz
    % Apply minimum image convention
    dxyz = mod(dxyz+Length/2,Length)-Length/2;
    B(:,j) = sqrt(sum(dxyz.^2,3));
end

% Generate graph of all chains in the system: nodes = chains

% Calculate adjacency matrix for sites. Do this by applying distance
% criteria to find adjacent sites and then converting to 1's and 0's:
C = ((B < 1.75) & (B > 0));
% Adjacency matrix for chains: N x N matrix such that nij = 1 if two chains
% are adjacent if any sites between them are adjacent and 0 otherwise:
L = StickSitesInSystem(:,3)+1; % Chain ID for each sticker site
S = sparse(L,1:numel(L),1); % sij = 1 if site is in chain, 0 otherwise
D = S*C*S' > 0; % Adjacency matrix except dii = 1
D(1:size(D,1)+1:end) = 0; % Adjacency matrix for chains
SystemChainsStickG = graph(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate graph of all chains in condensate: nodes = chains

% An edge in SystemChainsStickG means that at least one site on two chains
% is adjacent. Use something like union-find to find the largest sub-graph.
bincell = biconncomp(SystemChainsStickG, 'OutputForm', 'cell');
[s,d] = cellfun(@size,bincell);
out = max([s,d],[],'omitnan');
ind = d == out | s == out;
bincell_largest = bincell(ind);
ind_chain = cell2mat(bincell_largest); % Which nodes
% Generate largest sub-graph. Will be unique here.
% CondensateChainsStickG = subgraph(SystemChainsStickG,cell2mat(bincell_largest));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate graph of sticker-sticker interactions in condensate: nodes = sites

% Determine all sites in condensate
SitesInCondensate = find(ismember(A_all(:,3)+1,ind_chain));
% The above is with respect to system sites IDs.

% Modify graph such that the edges represent covalent bonding within each
% chain as well as non-bonded, sticker-sticker interactions. Discard
% non-bonded interactions that are not between two stickers. Treat all of
% the interactions depicted by the graph as being the same.

% Build adjacency matrix with edges as covalent bonds among sites. Note
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

% Separately, build adjacency matrix for sticker-sticker interactions in
% condensate and pad with zeros so that it is the same size as the
% adjacency matrix with respect to the covalent bonds, i.e., number of
% nodes = number of sites in the dense phase.

% Get adjacency matrix for stickers in condensate
StickSitesInSystem = find(ismember(A_all(:,2),Stickers)); % with respect to system sites IDs
StickSitesInCondensate = intersect(SitesInCondensate,StickSitesInSystem); % with respect to system sites IDs
StickSitesInDilute = find(~ismember(StickSitesInSystem,StickSitesInCondensate)); % with respect to sticker sites IDs
% Remove dilute-phase stickers from the sticker-sticker adjacency matrix
C(ismember(1:size(C,1),StickSitesInDilute),:) = [];
C(:,ismember(1:size(C,2),StickSitesInDilute)) = [];
% Pad adjacency matrix with zeros
stickerIDs = find(ismember(SitesInCondensate,StickSitesInCondensate)); % node IDs for stickers
C_stick = sparse(length(SitesInCondensate),length(SitesInCondensate));
C_stick(stickerIDs,stickerIDs) = C_stick(stickerIDs,stickerIDs) | C;

% Combine adjacency matrices
C_condensate = (C_bond | C_stick == 1);

% Generate graph for all sites in condensate, where edges represent
% non-bonded sticker-sticker interactions and chain connectivity
CondensateSitesG = graph(C_condensate);

% Remove self-loops
CondensateSitesG = simplify(CondensateSitesG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find site type for each site on single chain

siteType = A(1:Sites,2); % Will need to weight the sub-graphs