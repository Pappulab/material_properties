% Script to calculate viscoelastic properties of condensates and render
% single-chain sub-graphs

% Uses graph with chains as residues

% To run, this program needs generate_site_graph.m,
% calculate_single_chain_moduli.m, intersections.m, natsort.m, and
% natsortfiles.m, with the latter three programs from the MathWorks File
% Exchange. Takes as input data in LAMMPS trajectory format.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User input

N = 200; % Chains in system
Sites = 137; % Sites in each chain
Res = 11; % Number of unique residue types
Length = 120; % Length of cubic box on each side
SnapTot = 605; % Total number of snapshots in trajectory
SnapEq = 10; % Number of equilibrated snapshots to analyze at end of file (> 1)
Repl = 3; % Number of simulations (replicas) per temperature
b = 1; % Kuhn length
phi = 1; % Volume fraction (even though phi << 1, this doesn't matter)
xi = 1; % Friction coefficient of the background

% Relative energies

RelEl = zeros(Res,Res);

% These energies are specific to the A1-LCD system and are needed for
% weighting the sub-graphs. Make this an upper triangular array:

RelEl(1,1) = -3.15; % GG
RelEl(1,2) = -2.65; % GS
RelEl(1,3) = -2.35; % GT
RelEl(1,4) = -3.05; % GN
RelEl(1,5) = -2.75; % GQ
RelEl(1,6) = -3.15; % GF
RelEl(1,7) = -3.15; % GY
RelEl(1,8) = -3.15; % GR
RelEl(1,9) = -0.389; % GK
RelEl(1,10) = -3.15; % GD
RelEl(1,11) = -3.15; % GX

RelEl(2,2) = -2.65; % SS
RelEl(2,3) = -2.35; % ST
RelEl(2,4) = -3.05; % SN
RelEl(2,5) = -2.75; % SQ
RelEl(2,6) = -2.65; % SF
RelEl(2,7) = -2.65; % SY
RelEl(2,8) = -2.65; % SR
RelEl(2,9) = -0.389; % SK
RelEl(2,10) = -2.65; % SD
RelEl(2,11) = -2.65; % SX

RelEl(3,3) = -2.35; % TT
RelEl(3,4) = -2.35; % TN
RelEl(3,5) = -2.35; % TQ
RelEl(3,6) = -2.35; % TF
RelEl(3,7) = -2.35; % TY
RelEl(3,8) = -2.35; % TR
RelEl(3,9) = -0.389; % TK
RelEl(3,10) = -2.35; % TD
RelEl(3,11) = -2.35; % TX

RelEl(4,4) = -3.05; % NN
RelEl(4,5) = -2.75; % NQ
RelEl(4,6) = -3.05; % NF
RelEl(4,7) = -3.05; % NY
RelEl(4,8) = -3.05; % NR
RelEl(4,9) = -0.389; % NK
RelEl(4,10) = -3.05; % ND
RelEl(4,11) = -3.05; % NX

RelEl(5,5) = -2.75; % QQ
RelEl(5,6) = -2.75; % QF
RelEl(5,7) = -2.75; % QY
RelEl(5,8) = -2.75; % QR
RelEl(5,9) = -0.389; % QK
RelEl(5,10) = -2.75; % QD
RelEl(5,11) = -2.75; % QX

RelEl(6,6) = -13.4; % FF
RelEl(6,7) = -15.2; % FY
RelEl(6,8) = -11.0; % FR
RelEl(6,9) = -0.389; % FK
RelEl(6,10) = -2.95; % FD
RelEl(6,11) = -2.95; % FX

RelEl(7,7) = -19.3; % YY
RelEl(7,8) = -11.0; % YR
RelEl(7,9) = -0.389; % YK
RelEl(7,10) = -2.95; % YD
RelEl(7,11) = -2.95; % YX

RelEl(8,8) = -2.95; % RR
RelEl(8,9) = -0.389; % RK
RelEl(8,10) = -2.95; % RD
RelEl(8,11) = -2.95; % RX

RelEl(9,9) = -0.389; % KK
RelEl(9,10) = -0.389; % KD
RelEl(9,11) = -0.389; % KX

RelEl(10,10) = -2.95; % DD
RelEl(10,11) = -2.95; % DX

RelEl(11,11) = -2.95; % XX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analyze all .lammpstrj files in directory
datFiles = natsortfiles(dir('*.lammpstrj'));
numfiles = length(datFiles);

SitesInSystem = N*Sites;

avgGt = cell(numfiles,1);
avgStorage = cell(numfiles,1);
avgLoss = cell(numfiles,1);
avgViscosity = zeros(numfiles,1);
avgCompliance = zeros(numfiles,1);
tau_times = cell(numfiles,1);
std_tau_times = cell(numfiles,1);
avgX0 = zeros(numfiles,1);
avgY0 = zeros(numfiles,1);
stdGt = cell(numfiles,1);
stdStorage = cell(numfiles,1);
stdLoss = cell(numfiles,1);
stdViscosity = zeros(numfiles,1);
stdCompliance = zeros(numfiles,1);
stdX0 = zeros(numfiles,1);
stdY0 = zeros(numfiles,1);
repGraph = cell(numfiles,1);
time_vector = cell(numfiles,1);
freq_vector = cell(numfiles,1);

% Will need the following to take the continuous-time Fourier transform:
syms t w real

for k = 1:numfiles

    [~,baseFileNameNoExt,~] = fileparts(datFiles(k).name);

    % Truncate filename. Will need this later.
    TruncName = extractAfter(baseFileNameNoExt,'trj_');

    % Extract temperature. Will also need this later.
    filename = datFiles(k).name;
    subLoc = strfind(filename, '_');
    T = str2double(filename(subLoc(2) + 2 : subLoc(3) - 1));
    % The above assumes the filename is in the format: system_trj_TXX_
    % To do: make this more robust

    % Open and read the file
    fid = fopen(datFiles(k).name,'rt');
    chr = reshape(fread(fid,'*char'),1,[]);
    [~] = fclose(fid);
    cac = regexp(chr,'ITEM: TIMESTEP\n','split');
    len = size(cac,2);
    data = nan((len-1)*SitesInSystem,7);  % cac{1} is empty
    
    for jj = 2:len
        ccc = textscan(cac{jj},'%d%d%d%d%d%d%d','Headerlines',8,'CollectOutput',true);
        data(SitesInSystem*(jj-1)-SitesInSystem+1:SitesInSystem*(jj-1),:) = ccc{1};
    end

    % Skip beginning of trajectory
    data(1:(SnapTot-SnapEq)*SitesInSystem,:) = [];

    % Analyze each equilibrated snapshot
    out1 = cell(SnapEq,1);
    out2 = cell(SnapEq,1);
    out3 = cell(SnapEq,1);
    out4 = zeros(SnapEq,1);
    out5 = zeros(SnapEq,1);
    out6 = zeros(SnapEq,1);
    out7 = zeros(SnapEq,1);
    out8 = cell(SnapEq,1);
    out9 = cell(SnapEq,1);
    out10 = zeros(SnapEq,1);
    out11 = zeros(SnapEq,1);

    for n = 1:SnapEq
        A = data((SitesInSystem*n-SitesInSystem+1):SitesInSystem*n,:);

        % Generate relevant single-chain graphs
        [CondensateSitesG,ind_chain,siteType] ...
            = generate_site_graph(A,Length,Sites);

        % Calculate moduli
        [tvec,Gt,Freq,Storage,Loss,avg_Visc,std_Visc,avg_Comp,std_Comp,...
            tau_all,std_tau_all,X0,Y0,ChainsSubG] = ...
            calculate_single_chain_moduli(CondensateSitesG,ind_chain,...
            Sites,siteType,b,xi,T,phi,RelEl,t,w);

        % Relaxation modulus
        out1{n} = Gt;

        % Dynamic moduli
        out2{n} = Storage;
        out3{n} = Loss;

        % Viscosity
        out4(n) = avg_Visc;
        out5(n) = std_Visc;

        % Compliance
        out6(n) = avg_Comp;
        out7(n) = std_Comp;

        % Modes
        out8{n} = tau_all';
        out9{n} = std_tau_all';

        % Crossover
        out10(n) = X0;
        out11(n) = Y0;

    end

    % Time and frequency vectors should be the same for each snapshot at a
    % given temperature, but might be different at different temperatures.
    time_vector{k} = tvec;
    freq_vector{k} = Freq;

    % Average over all snapshots

    % Calculate average relaxation modulus per trajectory
    avgGt_snap = vertcat(mean(cell2mat(out1(:))));
    avgGt{k} = avgGt_snap;
    % Calculate standard deviation per trajectory
    stdGt_snap = vertcat(std(cell2mat(out1(:))));    
    stdGt{k} = stdGt_snap;

    % Calculate average storage modulus per trajectory
    avgStorage_snap = vertcat(mean(cell2mat(out2(:))));
    avgStorage{k} = avgStorage_snap;
    % Calculate standard deviation per trajectory
    stdStorage_snap = vertcat(std(cell2mat(out2(:))));    
    stdStorage{k} = stdStorage_snap;

    % Calculate average loss modulus per trajectory
    avgLoss_snap = vertcat(mean(cell2mat(out3(:))));
    avgLoss{k} = avgLoss_snap;
    % Calculate standard deviation per trajectory
    stdLoss_snap = vertcat(std(cell2mat(out3(:))));    
    stdLoss{k} = stdLoss_snap;

    % Calculate average viscosity per trajectory
    avgViscosity(k) = vertcat(mean(out4(:)));
    % Propagate error for average viscosity per trajectory
    varViscosity_snap = (1/SnapEq)^2*sum(reshape(out5(:).^2,SnapEq,[]));
    stdViscosity(k) = sqrt(varViscosity_snap);

    % Calculate average compliance per trajectory
    avgCompliance(k) = vertcat(mean(out6(:)));
    % Propagate error for average compliance per trajectory
    varCompliance_snap = (1/SnapEq)^2*sum(reshape(out7(:).^2,SnapEq,[]));
    stdCompliance(k) = sqrt(varCompliance_snap);

    % Calculate average of relaxation times per trajectory
    tau_times_snap = vertcat(mean(cell2mat(out8(:))));
    tau_times{k} = tau_times_snap;
    % Calculate standard deviation per trajectory
    std_tau_times_snap = vertcat(std(cell2mat(out9(:))));    
    std_tau_times{k} = std_tau_times_snap;

    % Calculate average crossover per trajectory
    avgX0(k) = vertcat(mean(out10(:)));
    avgY0(k) = vertcat(mean(out11(:)));
    % Calculate standard deviation per trajectory
    stdX0_snap = vertcat(std(out10(:)));
    stdY0_snap = vertcat(std(out11(:)));
    stdX0(k) = stdX0_snap;
    stdY0(k) = stdY0_snap;

    % Save representative chain graph using last snapshot in each trajectory
    repGraph{k} = ChainsSubG;
    deg = degree(repGraph{k}); % Degree

    % Generate and save unweighted graphs
    fig1 = figure('visible','off');
    plot(repGraph{k},'MarkerSize',5,'NodeCData',deg,'EdgeAlpha',1);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    set(gca,'XTick',[],'YTick',[])
    colorbar
    set(gca,'FontSize',24)
    str_append = '-single';
    str_saveas = sprintf('%s',TruncName,str_append);
    saveas(fig1,str_saveas,'png')

    % Save graph object as .mat file for further analysis
    newG = repGraph{k};
    save(str_saveas,'newG')

    close(fig1)

end

% Average over replicas for each temperature

% Moduli, viscosity, compliance, relaxation times, and crossover

relaxation_modulus = reshape(avgGt,Repl,[]);
% Compare replicas
minLength_t = length(tvec);
relaxation_modulus2 = reshape(mean(reshape(vertcat(relaxation_modulus{:}),Repl,[])),[],minLength_t);

storage_modulus = reshape(avgStorage,Repl,[]);
minLength_Freq = length(Freq); % Same for G"
storage_modulus2 = reshape(mean(reshape(vertcat(storage_modulus{:}),Repl,[])),[],minLength_Freq);

loss_modulus = reshape(avgLoss,Repl,[]);
loss_modulus2 = reshape(mean(reshape(vertcat(loss_modulus{:}),Repl,[])),[],minLength_Freq);

viscosity = mean(reshape(avgViscosity,Repl,[]));
compliance = mean(reshape(avgCompliance,Repl,[]));

relaxation_times = reshape(tau_times,Repl,[]);
minLength_times = Sites-1;
relaxation_times2 = reshape(mean(reshape(vertcat(relaxation_times{:}),Repl,[])),[],minLength_times);

crossover_x = mean(reshape(avgX0,Repl,[]));
crossover_y = mean(reshape(avgY0,Repl,[]));

% Propagate error for moduli

std_relaxation_modulus = reshape(stdGt,Repl,[]);
% Compare replicas and propagate error
var_relaxation_modulus2 = (1/Repl)^2*reshape(sum(reshape(vertcat(std_relaxation_modulus{:}).^2,Repl,[])),[],minLength_t);
std_relaxation_modulus2 = sqrt(var_relaxation_modulus2); % Standard deviation

std_storage_modulus = reshape(stdStorage,Repl,[]);
var_storage_modulus2 = (1/Repl)^2*reshape(sum(reshape(vertcat(std_storage_modulus{:}).^2,Repl,[])),[],minLength_Freq);
std_storage_modulus2 = sqrt(var_storage_modulus2); % Standard deviation

std_loss_modulus = reshape(stdLoss,Repl,[]);
var_loss_modulus2 = (1/Repl)^2*reshape(sum(reshape(vertcat(std_loss_modulus{:}).^2,Repl,[])),[],minLength_Freq);
std_loss_modulus2 = sqrt(var_loss_modulus2); % Standard deviation

% Propagate error for viscosity and compliance

varVisc = (1/Repl)^2*sum(reshape(stdViscosity.^2,Repl,[]));
Svisc = sqrt(varVisc); % Standard deviation of viscosity

varCompl = (1/Repl)^2*sum(reshape(stdCompliance.^2,Repl,[]));
Scompl = sqrt(varCompl); % Standard deviation of compliance

% Propagate error for relaxation times

std_relaxation_times = reshape(std_tau_times,Repl,[]);
var_relaxation_times2 = (1/Repl)^2*reshape(sum(reshape(vertcat(std_relaxation_times{:}),Repl,[])),[],minLength_times);
std_relaxation_times2 = sqrt(var_relaxation_times2); % Standard deviation

% Propagate error for crossover

varX0 = (1/Repl)^2*sum(reshape(stdX0.^2,Repl,[]));
SX0 = sqrt(varX0); % Standard deviation of X0

varY0 = (1/Repl)^2*sum(reshape(stdY0.^2,Repl,[]));
SY0 = sqrt(varY0); % Standard deviation of Y0

% Time and frequency vectors at each temperature

% This takes the mean for each replica, since the values can vary slightly.

time = reshape(time_vector,Repl,[]);
time2 = reshape(mean(reshape(vertcat(time{:}),Repl,[])),[],minLength_t);

frequency = reshape(freq_vector,Repl,[]);
frequency2 = reshape(mean(reshape(vertcat(frequency{:}),Repl,[])),[],minLength_Freq);

% Write to file

newC = cat(2,num2cell(time2',1),num2cell(relaxation_modulus2',1),...
    num2cell(std_relaxation_modulus2',1),num2cell(frequency2',1),...
    num2cell(storage_modulus2',1),num2cell(std_storage_modulus2',1),...
    num2cell(loss_modulus2',1),num2cell(std_loss_modulus2',1),...
    {viscosity',Svisc',compliance',Scompl',crossover_x',SX0',...
    crossover_y',SY0'})';
n2 = cellfun(@numel,newC);
n_max = max(n2);
n_vectors = numel(newC);
C_aug = cell(n_max,n_vectors);
for ii = 1:n_vectors
    C_aug(1:n2(ii),ii) = num2cell(newC{ii}(:));
end
writetable(cell2table(C_aug),'moduli_single_chain.csv','WriteVariableNames',0)

% Include modes

newC = cat(2,num2cell(relaxation_times2',1),num2cell(std_relaxation_times2',1))';
n2 = cellfun(@numel,newC);
n_max = max(max(n2));
n_vectors = numel(newC);
C_aug = cell(n_max,n_vectors);
for ii = 1:n_vectors
    C_aug(1:n2(ii),ii) = num2cell(newC{ii}(:));
end
writetable(cell2table(C_aug),'times_single_chain.csv','WriteVariableNames',0)

