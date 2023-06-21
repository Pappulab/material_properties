% Script to calculate viscoelastic properties of condensates

% Uses graph with chains as nodes; adjacencies are defined by stickers only

% To run, this program needs generate_sticker_chain_graph.m,
% calculate_overall_moduli.m, intersections.m, natsort.m, and
% natsortfiles.m, with the latter three programs from the MathWorks File
% Exchange. Takes as input data in LAMMPS trajectory format.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User input

N = 200; % Chains in system
Sites = 137; % Sites in each chain
Length = 120; % Length of cubic box on each side
SnapTot = 605; % Total number of snapshots in trajectory
SnapEq = 10; % Number of equilibrated snapshots to analyze at end of file (> 1)
Repl = 3; % Number of simulations (replicas) per temperature
b = 1; % Kuhn length
rho = 1; % Number of beads per unit volume
xi = 1; % Friction coefficient of the background
Stickers = [6,7]; % Sticker types

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
avgX0 = zeros(numfiles,1);
avgY0 = zeros(numfiles,1);
tau_times = cell(numfiles,1);
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
    out6 = cell(SnapEq,1);
    out7 = zeros(SnapEq,1);
    out8 = zeros(SnapEq,1);

    for n = 1:SnapEq
        A = data((SitesInSystem*n-SitesInSystem+1):SitesInSystem*n,:);

        % Generate relevant graphs
        CondensateChainsG = generate_sticker_chain_graph(A,Length,Stickers);

        % Calculate moduli
        [tvec,Gt,Freq,Storage,Loss,Visc,Comp,tau,X0,Y0]...
            = calculate_overall_moduli(CondensateChainsG,b,T,rho,xi,t,w);

        % Relaxation modulus
        out1{n} = Gt;

        % Dynamic moduli
        out2{n} = Storage;
        out3{n} = Loss;

        % Viscosity
        out4(n) = Visc;

        % Compliance
        out5(n) = Comp;

        % Modes
        out6{n} = tau;

        % Crossover
        out7(n) = X0;
        out8(n) = Y0;

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
    % Calculate standard deviation per trajectory
    stdViscosity_snap = vertcat(std(out4(:)));
    stdViscosity(k) = stdViscosity_snap;

    % Calculate average compliance per trajectory
    avgCompliance(k) = vertcat(mean(out5(:)));
    % Calculate standard deviation per trajectory
    stdCompliance_snap = vertcat(std(out5(:)));
    stdCompliance(k) = stdCompliance_snap;

    % Output modes (not average, since number of nodes in each graph might
    % not be the same)
    tau_times{k} = out6;

    % Calculate average crossover per trajectory
    avgX0(k) = vertcat(mean(out7(:)));
    avgY0(k) = vertcat(mean(out8(:)));
    % Calculate standard deviation per trajectory
    stdX0_snap = vertcat(std(out7(:)));
    stdY0_snap = vertcat(std(out8(:)));
    stdX0(k) = stdX0_snap;
    stdY0(k) = stdY0_snap;

    % Save representative graph using last snapshot in each trajectory
    repGraph{k} = CondensateChainsG;
    deg = degree(repGraph{k}); % Degree

    % Generate and save unweighted graphs
    fig1 = figure('visible','off');
    d = plot(repGraph{k},'MarkerSize',5,'NodeCData',deg,'EdgeAlpha',1);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    set(gca,'XTick',[],'YTick',[])
    colorbar
    set(gca,'FontSize',24)
    str_append = '-overall_A';
    str_saveas = sprintf('%s',TruncName,str_append);
    saveas(fig1,str_saveas,'png')

    % Save graph object as .mat file for further analysis
    newG = repGraph{k};
    save(str_saveas,'newG')

    close(fig1)

end

% Average over replicas for each temperature

% Moduli, viscosity, compliance, and crossover

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

crossover_x = mean(reshape(avgX0,Repl,[]));
crossover_y = mean(reshape(avgY0,Repl,[]));

% Propagate error for moduli

std_relaxation_modulus = reshape(stdGt,Repl,[]);
% Compare replicas and propagate error
std_relaxation_modulus2_sq = (1/Repl)^2*reshape(sum(reshape(vertcat(std_relaxation_modulus{:}).^2,Repl,[])),[],minLength_t);
std_relaxation_modulus2 = sqrt(std_relaxation_modulus2_sq); % Standard deviation

std_storage_modulus = reshape(stdStorage,Repl,[]);
std_storage_modulus2_sq = (1/Repl)^2*reshape(sum(reshape(vertcat(std_storage_modulus{:}).^2,Repl,[])),[],minLength_Freq);
std_storage_modulus2 = sqrt(std_storage_modulus2_sq); % Standard deviation

std_loss_modulus = reshape(stdLoss,Repl,[]);
std_loss_modulus2_sq = (1/Repl)^2*reshape(sum(reshape(vertcat(std_loss_modulus{:}).^2,Repl,[])),[],minLength_Freq);
std_loss_modulus2 = sqrt(std_loss_modulus2_sq); % Standard deviation

% Propagate error for viscosity and compliance

Sviscsq = (1/Repl)^2*sum(reshape(stdViscosity.^2,Repl,[]));
Svisc = sqrt(Sviscsq); % Standard deviation of viscosity

Scomplsq = (1/Repl)^2*sum(reshape(stdCompliance.^2,Repl,[]));
Scompl = sqrt(Scomplsq); % Standard deviation of compliance

% Propagate error for crossover

SX0sq = (1/Repl)^2*sum(reshape(stdX0.^2,Repl,[]));
SX0 = sqrt(SX0sq); % Standard deviation of X0

SY0sq = (1/Repl)^2*sum(reshape(stdY0.^2,Repl,[]));
SY0 = sqrt(SY0sq); % Standard deviation of Y0

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
writetable(cell2table(C_aug),'moduli_overall_A.csv','WriteVariableNames',0)

% Include modes

newC = cat(2,[tau_times{:}]')';
n2 = cellfun(@numel,newC);
n_max = max(max(n2));
n_vectors = numel(newC);
C_aug = cell(n_max,n_vectors);
for ii = 1:n_vectors
    C_aug(1:n2(ii),ii) = num2cell(newC{ii}(:));
end
writetable(cell2table(C_aug),'times_overall_A.csv','WriteVariableNames',0)
