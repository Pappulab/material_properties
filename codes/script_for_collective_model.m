% Script to calculate viscoelastic properties of condensates using
% Rouse-Zimm theory applied to a collective model as described in
% Alshareedah et al., Nature Physics 2024.

% Takes as input data in LAMMPS trajectory format. Multiple temperatures
% and replicates can be included, but only one type of sequence can be
% analyzed at a time. Note temperature is assumed to be in energy units.
% Also, note headers are not given in the output .csv files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User input

N = 200; % Chains in simulation
Sites = 137; % Sites in each chain
Length = 120; % Length of cubic box on each side
SnapTot = 605; % Total number of snapshots in trajectory
SnapEq = 10; % Number of snapshots to analyze at end of file (> 1)
Repl = 3; % Number of simulations (replicates) per temperature
phi = 1; % Volume fraction (set this to phi_dens / phi_sat as needed in post-processing)
b = 1; % Kuhn length
zeta = 1; % Friction coefficient of the background

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

% Take the Fourier transform symbolically
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

    % Analyze equilibrated snapshots
    out1 = cell(SnapEq,1);
    out2 = cell(SnapEq,1);
    out3 = cell(SnapEq,1);
    out4 = cell(SnapEq,1);
    out5 = cell(SnapEq,1);
    out6 = zeros(SnapEq,1);
    out7 = zeros(SnapEq,1);
    out8 = cell(SnapEq,1);
    out9 = zeros(SnapEq,1);
    out10 = zeros(SnapEq,1);

    for n = 1:SnapEq
        A = data((SitesInSystem*n-SitesInSystem+1):SitesInSystem*n,:);

        % Generate relevant graphs
        CondensateChainsG = generate_chain_graph(A,Length);

        % Calculate moduli
        [tvec,Gt,Freq,Storage,Loss,Visc,Comp,tau,X0,Y0]...
            = calculate_collective_moduli(CondensateChainsG,b,T,phi,zeta,t,w);

        % Time and frequency
        out1{n} = tvec;
        out2{n} = Freq;

        % Relaxation modulus
        out3{n} = Gt;

        % Dynamic moduli
        out4{n} = Storage;
        out5{n} = Loss;

        % Viscosity
        out6(n) = Visc;

        % Compliance
        out7(n) = Comp;

        % Relaxation times
        out8{n} = tau;

        % Crossover
        out9(n) = X0;
        out10(n) = Y0;

    end

    % Average over all snapshots per replicate

    % Relaxation modulus
    [time_vector{k},avgGt{k}] = calculate_mean_of_cell_array_and_concat(out1,out3,SnapEq);
    % Standard deviation
    stdGt{k} = calculate_std_of_cell_array_and_concat(out1,out3,SnapEq);

    % Storage modulus
    [freq_vector{k},avgStorage{k}] = calculate_mean_of_cell_array_and_concat(out2,out4,SnapEq);
    % Standard deviation
    stdStorage{k} = calculate_std_of_cell_array_and_concat(out2,out4,SnapEq);

    % Loss modulus
    [~,avgLoss{k}] = calculate_mean_of_cell_array_and_concat(out2,out5,SnapEq);
    % Standard deviation
    stdLoss{k} = calculate_std_of_cell_array_and_concat(out2,out5,SnapEq);  

    % Viscosity
    avgViscosity(k) = calculate_mean_and_concat(out6);
    stdViscosity(k) = calculate_std_and_concat(out6);

    % Compliance
    avgCompliance(k) = calculate_mean_and_concat(out7);
    stdCompliance(k) = calculate_std_and_concat(out7);

    % Output relaxation times (not average, since number of nodes in each
    % graph using the collective model might not be the same)
    tau_times{k} = out8;

    % Crossover
    avgX0(k) = calculate_mean_and_concat(out9);
    avgY0(k) = calculate_mean_and_concat(out10);
    % Standard deviation
    stdX0(k) = calculate_std_and_concat(out9);
    stdY0(k) = calculate_std_and_concat(out10);

    % Save representative graph using last snapshot in each trajectory
    repGraph{k} = CondensateChainsG;
    deg = degree(repGraph{k}); % Degree

    % Generate and save graphs
    fig1 = figure('visible','off');
    plot(repGraph{k},'MarkerSize',5,'NodeCData',deg,'EdgeAlpha',1);
    set(gca,'XTickLabel',[],'YTickLabel',[])
    set(gca,'XTick',[],'YTick',[])
    axis off
    colorbar
    C = parula(max(deg)-min(deg)+1); % Want discrete values for the degree
    colormap(C)
    set(gca,'FontSize',24)
    str_append = '-collective';
    str_saveas = sprintf('%s',TruncName,str_append);
    saveas(fig1,str_saveas,'png')

    % Save graph object as .mat file for further analysis
    newG = repGraph{k};
    save(str_saveas,'newG')

    close(fig1)

end

% Average over replicates for each temperature

% Relaxation modulus
[time,relaxation_modulus] = calculate_mean_over_replicates(time_vector,avgGt,Repl,numfiles);
% Standard deviation
std_relaxation_modulus = propagate_std_vector_over_replicates(time_vector,stdGt,Repl,numfiles);

% Dynamical moduli
[frequency,storage_modulus] = calculate_mean_over_replicates(freq_vector,avgStorage,Repl,numfiles);
[~,loss_modulus] = calculate_mean_over_replicates(freq_vector,avgLoss,Repl,numfiles);
% Standard deviation
std_storage_modulus = propagate_std_vector_over_replicates(freq_vector,stdStorage,Repl,numfiles);
std_loss_modulus = propagate_std_vector_over_replicates(freq_vector,stdLoss,Repl,numfiles);

% Viscosity
viscosity = mean(reshape(avgViscosity,Repl,[]));
% Standard devation
std_visc = propagate_std_value_over_replicates(stdViscosity,Repl);

% Compliance
compliance = mean(reshape(avgCompliance,Repl,[]));
% Standard deviation
std_compl = propagate_std_value_over_replicates(stdCompliance,Repl);

% Crossover
crossover_x = mean(reshape(avgX0,Repl,[]));
crossover_y = mean(reshape(avgY0,Repl,[]));

% Standard deviation in crossover
std_X0 = propagate_std_value_over_replicates(stdX0,Repl);
std_Y0 = propagate_std_value_over_replicates(stdY0,Repl);

% Write to files
write_rheology_data_to_csv(time, relaxation_modulus, std_relaxation_modulus, ...
    frequency, storage_modulus, std_storage_modulus, loss_modulus, ...
    std_loss_modulus, viscosity, std_visc, compliance, std_compl, crossover_x, ...
    std_X0, crossover_y, std_Y0);
write_times_to_csv(tau_times);

function CondensateChainsG = generate_chain_graph(A,Length)
    
    % Makes graph of all chains in condensate using physical contact between
    % residues to define chain adjacencies. (This could be made a lot more
    % efficient.)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find which chains are in the condensate by first making graph of all
    % chains in the system and then determining largest connected network. By
    % system, we mean dilute + dense phase.
    
    % 1. Calculate all possible pairwise distances between sites. Scales as O(N^2).
    
    SitesInSystem = length(A); % Number of sites in system
    
    Ai = reshape(A(:,4:6),[],1,3);
    B = zeros(SitesInSystem,SitesInSystem);
    for j = 1:SitesInSystem
        Aj = reshape(A(j,4:6),1,[],3);
        dxyz = Ai-Aj; % Differences dx, dy, dz
        % Apply minimum image convention
        dxyz = mod(dxyz+Length/2,Length)-Length/2; % Coordinates initially are [0,Length]
        B(:,j) = sqrt(sum(dxyz.^2,3));
    end
    
    % 2. Generate graph of all chains in the system: nodes = chains
    
    % Calculate adjacency matrix for sites. Do this by applying distance
    % criteria to find adjacent sites and then converting to 1's and 0's:
    C = ((B < 1.75) & (B > 0));
    % Adjacency matrix for chains: N x N matrix s.t. nij = 1 if two chains are
    % adjacent if any sites between them are adjacent and 0 otherwise:
    L = A(:,3)+1; % Chain ID for each site
    S = sparse(L,1:numel(L),1); % sij = 1 if site j is in chain i, 0 otherwise
    D = S*C*S' > 0; % Adjacency matrix except dii = 1
    D(1:size(D,1)+1:end) = 0; % Adjacency matrix for chains removing self-loops
    SystemChainsAllG = graph(D);
    
    % 3. Generate sub-graph of all chains in condensate: nodes = chains
    
    % An edge in SystemChainsAllG means that at least one site on two
    % chains is adjacent. Find the largest sub-graph.
    bincell = biconncomp(SystemChainsAllG, 'OutputForm', 'cell');
    [s,d] = cellfun(@size,bincell);
    out = max([s,d],[],'omitnan');
    ind = d == out | s == out;
    bincell_largest = bincell(ind);
    % Generate largest sub-graph. Will be unique here.
    CondensateChainsG = subgraph(SystemChainsAllG,cell2mat(bincell_largest));

end

function [tvec,Gt,Freq,Storage,Loss,Visc,Comp,tau,X0,Y0] = ...
    calculate_collective_moduli(CondensateChainsG,b,T,phi,zeta,t,w)
    
    % Calculates viscoelastic properties of condensate according to
    % graph-theoretic formulation of Rouse-Zimm model
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Number of nodes in graph
    NumChains = numnodes(CondensateChainsG);
    
    % Graph Laplacian aka Zimm matrix
    Zimm = laplacian(CondensateChainsG);
    
    % Find eigenvalues of Zimm matrix
    lambda = eig(Zimm);
    
    % Calculate relaxation times, treating simulation T in energy units
    tau = zeta*b^2*(6*T*lambda(2:end)).^-1; % First eigenvalue is zero
    
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
    % See equations in Rouse 1953. Effectively, this is the same as doing
    % the calculation symbolically.)
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
    Freq = logspace(min_order_w-3,max_order_w,1000); % All frequencies and then some
    Storage = double(Gw1(Freq));
    Loss = double(Gw2(Freq));
    
    % Calculate crossover
    [X0,Y0] = intersections(Freq,Storage,Freq,Loss,'robust');

end

% Bookkeeping

function result_mean_value = calculate_mean_and_concat(data)
    % Calculate the mean and concatenate the values
    result_mean_value = vertcat(mean(data(:)));
end

function result_std_value = calculate_std_and_concat(data)
    % Calculate the std and concatenate the values
    result_std_value = vertcat(std(data(:)));
end

function [x_out,y_out] = calculate_mean_of_cell_array_and_concat(x_in,y_in,SnapEq)
    % This assumes a little variation within a replicate

    % Generate common domain using unique x values from all vectors
    commonDomain = unique(vertcat(x_in{:}));
    
    % Preallocate variables
    sumYValues = zeros(numel(commonDomain), 1);
    validCounts = zeros(numel(commonDomain), 1);
    
    % Loop through each set of vectors
    for n = 1:SnapEq
        % Interpolate y values to match the common domain
        interpolatedYValues = interp1(x_in{n}, y_in{n}, commonDomain, 'linear');
        
        % Check for NaN values in the interpolated y values
        nanIndices = isnan(interpolatedYValues);
        
        % Exclude interpolated values that are NaN from the sum and count
        validIndices = ~nanIndices;
        
        % Add the valid interpolated y values to the sum
        sumYValues(validIndices) = sumYValues(validIndices) + interpolatedYValues(validIndices);
        
        % Increment the count of valid values
        validCounts(validIndices) = validCounts(validIndices) + 1;
    end
    
    % Calculate the average y values, accounting for cases with NaN
    y_out = sumYValues ./ max(validCounts, 1);

    x_out = commonDomain;
end

function result_std_vector = calculate_std_of_cell_array_and_concat(x_in, y_in, SnapEq)
    % Generate common domain using unique x values from all vectors
    commonDomain = unique(vertcat(x_in{:}));
    
    % Preallocate array to store all interpolated y values
    interpolatedYValues = zeros(numel(commonDomain), SnapEq);
    
    % Loop through each set of vectors
    for n = 1:SnapEq
        % Interpolate y values to match the common domain
        interpolatedYValues(:, n) = interp1(x_in{n}, y_in{n}, commonDomain, 'linear');
    end
    
    % Compute the standard deviation for each element across all sets of vectors
    result_std_vector = std(interpolatedYValues, 0, 2, 'omitnan'); % Compute along the second dimension with normalization by N-1
end

function [x_out,y_out] = calculate_mean_over_replicates(x_in,y_in,Repl,numfiles)
    % This assumes a little variation within replicates

    % Generate common domain using unique x values from all vectors
    allX = vertcat(x_in{:});
    commonDomain = unique(allX);
    
    % Initialize cell arrays to store interpolated data for each replicate and temperature
    interpolatedYValues = cell(Repl * numfiles/Repl, 1);
    
    % Interpolate each replicate onto the common domain for all temperatures
    for temp = 1:numfiles/Repl
        for i = 1:Repl
            % Index to access the current replicate's data
            index = (temp - 1) * Repl + i;
            
            % Interpolate data for the current replicate onto the common domain
            interpolatedYValues{index} = interp1(x_in{index}, y_in{index}, commonDomain, 'linear');
        end
    end

    % Preallocate cell arrays to store mean values for each temperature
    y_out = cell(1, numfiles/Repl);
    
    % Calculate the mean for each temperature and reshape
    for temp = 1:numfiles/Repl
        % Compute the mean for the replicates at the current temperature
        tempIndices = (temp - 1) * Repl + (1:Repl);

        % Concatenate the subsets into a 3D array
        subset_array = cat(3, interpolatedYValues{tempIndices});
        
        % Find NaN values in the array
        nan_indices = isnan(subset_array);
        
        % Replace NaN values with zero
        subset_array(nan_indices) = 0;
        
        % Count the number of non-NaN values along the third dimension
        non_nan_counts = sum(~nan_indices, 3);
        
        % Calculate the sum of non-NaN values along the third dimension
        sum_without_nans = sum(subset_array, 3);

        % Calculate the mean while ignoring NaNs
        mean_interpolated_data = sum_without_nans ./ non_nan_counts;

        % Reshape the mean data to match the length of the replicates
        length_replicate = size(interpolatedYValues{tempIndices(1)}, 1);
        y_out{temp} = reshape(mean_interpolated_data, [], length_replicate);
    end
    
    % Assign common domain for all temperatures
    x_out = commonDomain;

end

function std_out_propagated = propagate_std_vector_over_replicates(x_in,y_in,Repl,numfiles)
    % Generate common domain using unique x values from all vectors
    allX = vertcat(x_in{:});
    commonDomain = unique(allX);
    
    % Initialize cell array to store interpolated standard deviation data for each replicate and temperature
    interpolated_std_replicates = cell(Repl * numfiles/Repl, 1);
    
    % Interpolate each replicate's standard deviation onto the common domain for all temperatures
    for temp = 1:numfiles/Repl
        for i = 1:Repl
            % Perform interpolation
            interpolated_data = interp1(x_in{(temp-1)*Repl + i}, y_in{(temp-1)*Repl + i}, commonDomain, 'linear');

            % Handle NaN values
            nan_indices = isnan(interpolated_data);
            interpolated_data(nan_indices) = 0;

            % Perform linear interpolation on non-NaN values
            interpolated_std_replicates{(temp-1)*Repl + i} = interpolated_data;
        end
    end    
    transposed_std_replicates = cellfun(@transpose, interpolated_std_replicates, 'UniformOutput', false);

    % Error propagation for standard deviation
    var_result = (1/Repl)^2 * reshape(sum(reshape(vertcat(transposed_std_replicates{:}).^2, Repl, [])), [], size(interpolated_std_replicates{1}, 1));

    % Calculate the standard deviation
    std_out_propagated = sqrt(var_result);
end

function result_std_propaged_value = propagate_std_value_over_replicates(input_data, Repl)
    % Reshape the input data
    reshaped_data = reshape(input_data, Repl, []);

    % Error propagation for standard deviation
    var_result = (1/Repl)^2 * sum(reshape(reshaped_data.^2, Repl, []));

    % Calculate the standard deviation
    result_std_propaged_value = sqrt(var_result);
end

function write_rheology_data_to_csv(time, relaxation_modulus, std_relaxation_modulus, ...
    frequency, storage_modulus, std_storage_modulus, loss_modulus, ...
    std_loss_modulus, viscosity, std_visc, compliance, std_compl, crossover_x, std_X0, crossover_y, std_Y0)

    % Outputs the rheology data - messy, but it works!

    % Combine data into cell arrays
    transposed_std_relaxation_modulus = std_relaxation_modulus';
    transposed_std_storage_modulus = std_storage_modulus';
    transposed_std_loss_modulus = std_loss_modulus';    
    % Convert each column to a cell and concatenate them vertically
    std_Gt = cell(size(transposed_std_relaxation_modulus, 2), 1);
    std_Storage = cell(size(transposed_std_storage_modulus, 2), 1);
    std_Loss = cell(size(transposed_std_loss_modulus, 2), 1);
    for i = 1:size(transposed_std_relaxation_modulus, 2)
        std_Gt{i} = transposed_std_relaxation_modulus(:, i);
    end
    for i = 1:size(transposed_std_storage_modulus, 2)
        std_Storage{i} = transposed_std_storage_modulus(:, i);
        std_Loss{i} = transposed_std_loss_modulus(:, i);
    end

    newC = [{time}; cellfun(@(x) x(:), relaxation_modulus, 'UniformOutput', false).'; ...
        std_Gt; ...
        {frequency}; cellfun(@(x) x(:), storage_modulus, 'UniformOutput', false).'; ...
        std_Storage; ...
        cellfun(@(x) x(:), loss_modulus, 'UniformOutput', false).'; ...
        std_Loss; ...
        {viscosity';std_visc';compliance';std_compl';crossover_x';std_X0';crossover_y';std_Y0'}];

    % Calculate dimensions
    n2 = cellfun(@numel, newC);
    n_max = max(n2);
    n_vectors = numel(newC);

    % Augment the cell array to handle varying lengths
    C_aug = cell(n_max, n_vectors);
    for ii = 1:n_vectors
        C_aug(1:n2(ii), ii) = num2cell(newC{ii}(:));
    end

    % Write to .csv file
    writetable(cell2table(C_aug), 'moduli_collective.csv', 'WriteVariableNames', 0);
end

function write_times_to_csv(tau_times)
    % Combine times data into a cell array
    newC = cat(2, [tau_times{:}]')';

    % Calculate dimensions
    n2 = cellfun(@numel, newC);
    n_max = max(max(n2));
    n_vectors = numel(newC);

    % Augment the cell array to handle varying lengths
    C_aug = cell(n_max, n_vectors);
    for ii = 1:n_vectors
        C_aug(1:n2(ii), ii) = num2cell(newC{ii}(:));
    end

    % Write to .csv file
    writetable(cell2table(C_aug), 'times_collective.csv', 'WriteVariableNames', 0);
end

% From “Natural-Order Filename Sort” v3.4.5 on MathWorks File Exchange

function [B,ndx,dbg] = natsort(A,rgx,varargin)
    % Natural-order / alphanumeric sort the elements of a text array.
    %
    % (c) 2012-2023 Stephen Cobeldick
    %
    % Sorts text by character code and by number value. By default matches
    % integer substrings and performs a case-insensitive ascending sort.
    % Options to select the number format, sort order, case sensitivity, etc.
    %
    %%% Example:
    % >> A = ["x2", "x10", "x1"];
    % >> natsort(A)
    % ans =   "x1"  "x2"  "x10"
    %
    %%% Syntax:
    %  B = natsort(A)
    %  B = natsort(A,rgx)
    %  B = natsort(A,rgx,<options>)
    % [B,ndx,dbg] = natsort(A,...)
    %
    % To sort any file-names or folder-names use NATSORTFILES (File Exchange 47434)
    % To sort the rows of a string/cell/table use NATSORTROWS (File Exchange 47433)
    % To sort string/cells using custom sequences use ARBSORT (File Exchange 132263)
    %
    % Number Format %%
    %
    % The **default regular expression '\d+' matches consecutive digit
    % characters, i.e. integer numbers. Specifying the optional regular
    % expression allows the numbers to include a +/- sign, decimal point,
    % decimal fraction digits, exponent E-notation, character quantifiers,
    % or lookarounds. For information on defining regular expressions:
    % <http://www.mathworks.com/help/matlab/matlab_prog/regular-expressions.html>
    % For example, to match leading/trailing whitespace prepend/append '\s*'.
    %
    % The number substrings are parsed by SSCANF into numeric values, using
    % either the **default format '%f' or the user-supplied format specifier.
    % Both decimal comma and decimal point are accepted in number substrings.
    %
    % This table shows examples of some regular expression patterns for common
    % notations and ways of writing numbers, together with suitable SSCANF formats:
    %
    % Regular       | Number Substring | Number Substring              | SSCANF
    % Expression:   | Match Examples:  | Match Description:            | Format Specifier:
    % ==============|==================|===============================|==================
    % **        \d+ | 0, 123, 4, 56789 | unsigned integer              | %f  %i  %u  %lu
    % --------------|------------------|-------------------------------|------------------
    %      [-+]?\d+ | +1, 23, -45, 678 | integer with optional +/- sign| %f  %i  %d  %ld
    % --------------|------------------|-------------------------------|------------------
    %     \d+\.?\d* | 012, 3.45, 678.9 | integer or decimal            | %f
    % (\d+|Inf|NaN) | 123, 4, NaN, Inf | integer, Inf, or NaN          | %f
    %  \d+\.\d+E\d+ | 0.123e4, 5.67e08 | exponential notation          | %f
    % --------------|------------------|-------------------------------|------------------
    %  0X[0-9A-F]+  | 0X0, 0X3E7, 0XFF | hexadecimal notation & prefix | %x  %i
    %    [0-9A-F]+  |   0,   3E7,   FF | hexadecimal notation          | %x
    % --------------|------------------|-------------------------------|------------------
    %  0[0-7]+      | 012, 03456, 0700 | octal notation & prefix       | %o  %i
    %   [0-7]+      |  12,  3456,  700 | octal notation                | %o
    % --------------|------------------|-------------------------------|------------------
    %  0B[01]+      | 0B1, 0B101, 0B10 | binary notation & prefix      | %b   (not SSCANF)
    %    [01]+      |   1,   101,   10 | binary notation               | %b   (not SSCANF)
    % --------------|------------------|-------------------------------|------------------
    %
    % Debugging Output Array %%
    %
    % The third output is a cell array <dbg>, for checking how the numbers
    % were matched by the regular expression <rgx> and converted to numeric
    % by the SSCANF format. The rows of <dbg> are linearly indexed from
    % the first input argument <A>.
    %
    % >> [~,~,dbg] = natsort(A)
    % dbg =
    %    'x'    [ 2]
    %    'x'    [10]
    %    'x'    [ 1]
    %
    % Examples %%
    %
    %%% Multiple integers (e.g. release version numbers):
    % >> Aa = {'v10.6', 'v9.10', 'v9.5', 'v10.10', 'v9.10.20', 'v9.10.8'};
    % >> sort(Aa) % for comparison.
    % ans =   'v10.10'  'v10.6'  'v9.10'  'v9.10.20'  'v9.10.8'  'v9.5'
    % >> natsort(Aa)
    % ans =   'v9.5'  'v9.10'  'v9.10.8'  'v9.10.20'  'v10.6'  'v10.10'
    %
    %%% Integer, decimal, NaN, or Inf numbers, possibly with +/- signs:
    % >> Ab = {'test+NaN', 'test11.5', 'test-1.4', 'test', 'test-Inf', 'test+0.3'};
    % >> sort(Ab) % for comparison.
    % ans =   'test' 'test+0.3' 'test+NaN' 'test-1.4' 'test-Inf' 'test11.5'
    % >> natsort(Ab, '[-+]?(NaN|Inf|\d+\.?\d*)')
    % ans =   'test' 'test-Inf' 'test-1.4' 'test+0.3' 'test11.5' 'test+NaN'
    %
    %%% Integer or decimal numbers, possibly with an exponent:
    % >> Ac = {'0.56e007', '', '43E-2', '10000', '9.8'};
    % >> sort(Ac) % for comparison.
    % ans =   ''  '0.56e007'  '10000'  '43E-2'  '9.8'
    % >> natsort(Ac, '\d+\.?\d*(E[-+]?\d+)?')
    % ans =   ''  '43E-2'  '9.8'  '10000'  '0.56e007'
    %
    %%% Hexadecimal numbers (with '0X' prefix):
    % >> Ad = {'a0X7C4z', 'a0X5z', 'a0X18z', 'a0XFz'};
    % >> sort(Ad) % for comparison.
    % ans =   'a0X18z'  'a0X5z'  'a0X7C4z'  'a0XFz'
    % >> natsort(Ad, '0X[0-9A-F]+', '%i')
    % ans =   'a0X5z'  'a0XFz'  'a0X18z'  'a0X7C4z'
    %
    %%% Binary numbers:
    % >> Ae = {'a11111000100z', 'a101z', 'a000000000011000z', 'a1111z'};
    % >> sort(Ae) % for comparison.
    % ans =   'a000000000011000z'  'a101z'  'a11111000100z'  'a1111z'
    % >> natsort(Ae, '[01]+', '%b')
    % ans =   'a101z'  'a1111z'  'a000000000011000z'  'a11111000100z'
    %
    %%% Case sensitivity:
    % >> Af = {'a2', 'A20', 'A1', 'a10', 'A2', 'a1'};
    % >> natsort(Af, [], 'ignorecase') % default
    % ans =   'A1'  'a1'  'a2'  'A2'  'a10'  'A20'
    % >> natsort(Af, [], 'matchcase')
    % ans =   'A1'  'A2'  'A20'  'a1'  'a2'  'a10'
    %
    %%% Sort order:
    % >> Ag = {'2', 'a', '', '3', 'B', '1'};
    % >> natsort(Ag, [], 'ascend') % default
    % ans =   ''   '1'  '2'  '3'  'a'  'B'
    % >> natsort(Ag, [], 'descend')
    % ans =   'B'  'a'  '3'  '2'  '1'  ''
    % >> natsort(Ag, [], 'num<char') % default
    % ans =   ''   '1'  '2'  '3'  'a'  'B'
    % >> natsort(Ag, [], 'char<num')
    % ans =   ''   'a'  'B'  '1'  '2'  '3'
    %
    %%% UINT64 numbers (with full precision):
    % >> natsort({'a18446744073709551615z', 'a18446744073709551614z'}, [], '%lu')
    % ans =       'a18446744073709551614z'  'a18446744073709551615z'
    %
    % Input and Output Arguments %%
    %
    %%% Inputs (**=default):
    % A   = Array to be sorted. Can be a string array, or a cell array of
    %       character row vectors, or a categorical array, or a datetime array,
    %       or any other array type which can be converted by CELLSTR.
    % rgx = Optional regular expression to match number substrings.
    %     = [] uses the default regular expression '\d+'** to match integers.
    % <options> can be entered in any order, as many as required:
    %     = Sort direction: 'descend'/'ascend'**
    %     = Character case handling: 'matchcase'/'ignorecase'**
    %     = Character/number order: 'char<num'/'num<char'**
    %     = NaN/number order: 'NaN<num'/'num<NaN'**
    %     = SSCANF conversion format: e.g. '%x', '%li', '%b', '%f'**, etc.
    %     = Function handle of a function that sorts text. It must accept one
    %       input, which is a cell array of char vectors (the text array to
    %       be sorted). It must return as its 2nd output the sort indices.
    %
    %%% Outputs:
    % B   = Array <A> sorted into natural sort order.     The same size as <A>.
    % ndx = NumericArray, generally such that B = A(ndx). The same size as <A>.
    % dbg = CellArray of the parsed characters and number values. Each row
    %       corresponds to one input element of <A>, in linear-index order.
    %
    % See also SORT NATSORT_TEST NATSORTFILES NATSORTROWS ARBSORT
    % IREGEXP REGEXP COMPOSE STRING STRINGS CATEGORICAL CELLSTR SSCANF
    % Input Wrangling %%
    %
    fnh = @(c)cellfun('isclass',c,'char') & cellfun('size',c,1)<2 & cellfun('ndims',c)<3;
    %
    if iscell(A)
	    assert(all(fnh(A(:))),...
		    'SC:natsort:A:CellInvalidContent',...
		    'First input <A> cell array must contain only character row vectors.')
	    C = A(:);
    elseif ischar(A) % Convert char matrix:
	    assert(ndims(A)<3,...
		    'SC:natsort:A:CharNotMatrix',...
		    'First input <A> if character class must be a matrix.') %#ok<ISMAT>
	    C = num2cell(A,2);
    else % Convert string, categorical, datetime, enumeration, etc.:
	    C = cellstr(A(:));
    end
    %
    chk = '(match|ignore)(case|dia)|(de|a)scend(ing)?|(char|nan|num)[<>](char|nan|num)|%[a-z]+';
    %
    if nargin<2 || isnumeric(rgx)&&isequal(rgx,[])
	    rgx = '\d+';
    elseif ischar(rgx)
	    assert(ndims(rgx)<3 && size(rgx,1)==1,...
		    'SC:natsort:rgx:NotCharVector',...
		    'Second input <rgx> character row vector must have size 1xN.') %#ok<ISMAT>
	    nsChkRgx(rgx,chk)
    else
	    rgx = ns1s2c(rgx);
	    assert(ischar(rgx),...
		    'SC:natsort:rgx:InvalidType',...
		    'Second input <rgx> must be a character row vector or a string scalar.')
	    nsChkRgx(rgx,chk)
    end
    %
    varargin = cellfun(@ns1s2c, varargin, 'UniformOutput',false);
    ixv = fnh(varargin); % char
    txt = varargin(ixv); % char
    xtx = varargin(~ixv); % not
    %
    % Sort direction:
    tdd = strcmpi(txt,'descend');
    tdx = strcmpi(txt,'ascend')|tdd;
    % Character case:
    tcm = strcmpi(txt,'matchcase');
    tcx = strcmpi(txt,'ignorecase')|tcm;
    % Char/num order:
    ttn = strcmpi(txt,'num>char')|strcmpi(txt,'char<num');
    ttx = strcmpi(txt,'num<char')|strcmpi(txt,'char>num')|ttn;
    % NaN/num order:
    ton = strcmpi(txt,'num>NaN')|strcmpi(txt,'NaN<num');
    tox = strcmpi(txt,'num<NaN')|strcmpi(txt,'NaN>num')|ton;
    % SSCANF format:
    tsf = ~cellfun('isempty',regexp(txt,'^%([bdiuoxfeg]|l[diuox])$'));
    %
    nsAssert(txt, tdx, 'SortDirection', 'sort direction')
    nsAssert(txt, tcx,  'CaseMatching', 'case sensitivity')
    nsAssert(txt, ttx,  'CharNumOrder', 'number-character order')
    nsAssert(txt, tox,   'NanNumOrder', 'number-NaN order')
    nsAssert(txt, tsf,  'sscanfFormat', 'SSCANF format')
    %
    ixx = tdx|tcx|ttx|tox|tsf;
    if ~all(ixx)
	    error('SC:natsort:InvalidOptions',...
		    ['Invalid options provided. Check the help and option spelling!',...
		    '\nThe provided options:%s'],sprintf(' "%s"',txt{~ixx}))
    end
    %
    % SSCANF format:
    if any(tsf)
	    fmt = txt{tsf};
    else
	    fmt = '%f';
    end
    %
    xfh = cellfun('isclass',xtx,'function_handle');
    assert(nnz(xfh)<2,...
	    'SC:natsort:FunctionHandle:Overspecified',...
	    'The function handle option may only be specified once.')
    assert(all(xfh),...
	    'SC:natsort:InvalidOptions',...
	    'Optional arguments must be character row vectors, string scalars, or function handles.')
    if any(xfh)
	    txfh = xtx{xfh};
    end
    %
    % Identify and Convert Numbers %%
    %
    [nbr,spl] = regexpi(C(:), rgx, 'match','split', txt{tcx});
    %
    if numel(nbr)
	    V = [nbr{:}];
	    if strcmp(fmt,'%b')
		    V = regexprep(V,'^0[Bb]','');
		    vec = cellfun(@(s)pow2(numel(s)-1:-1:0)*sscanf(s,'%1d'),V);
	    else
		    vec = sscanf(strrep(sprintf(' %s','0',V{:}),',','.'),fmt);
		    vec = vec(2:end); % SSCANF wrong data class bug (R2009b and R2010b)
	    end
	    assert(numel(vec)==numel(V),...
		    'SC:natsort:sscanf:TooManyValues',...
		    'The "%s" format must return one value for each input number.',fmt)
    else
	    vec = [];
    end
    %
    % Allocate Data %%
    %
    % Determine lengths:
    nmx = numel(C);
    lnn = cellfun('length',nbr);
    lns = cellfun('length',spl);
    mxs = max(lns);
    %
    % Allocate data:
    idn = permute(bsxfun(@le,1:mxs,lnn),[2,1]); % TRANSPOSE lost class bug (R2013b)
    ids = permute(bsxfun(@le,1:mxs,lns),[2,1]); % TRANSPOSE lost class bug (R2013b)
    arn = zeros(mxs,nmx,class(vec));
    ars =  cell(mxs,nmx);
    ars(:) = {''};
    ars(ids) = [spl{:}];
    arn(idn) = vec;
    %
    % Debugging Array %%
    %
    if nargout>2
	    dbg = cell(nmx,0);
	    for k = 1:nmx
		    V = spl{k};
		    V(2,:) = [num2cell(arn(idn(:,k),k));{[]}];
		    V(cellfun('isempty',V)) = [];
		    dbg(k,1:numel(V)) = V;
	    end
    end
    %
    % Sort Matrices %%
    %
    if ~any(tcm) % ignorecase
	    ars = lower(ars);
    end
    %
    if any(ttn) % char<num
	    % Determine max character code:
	    mxc = 'X';
	    tmp = warning('off','all');
	    mxc(1) = Inf;
	    warning(tmp)
	    mxc(mxc==0) = 255; % Octave
	    % Append max character code to the split text:
	    %ars(idn) = strcat(ars(idn),mxc); % slower than loop
	    for ii = reshape(find(idn),1,[])
		    ars{ii}(1,end+1) = mxc;
	    end
    end
    %
    idn(isnan(arn)) = ~any(ton); % NaN<num
    %
    if any(xfh) % external text-sorting function
	    [~,ndx] = txfh(ars(mxs,:));
	    for ii = mxs-1:-1:1
		    [~,idx] = sort(arn(ii,ndx),txt{tdx});
		    ndx = ndx(idx);
		    [~,idx] = sort(idn(ii,ndx),txt{tdx});
		    ndx = ndx(idx);
		    [~,idx] = txfh(ars(ii,ndx));
		    ndx = ndx(idx);
	    end
    elseif any(tdd)
	    [~,ndx] = sort(nsGroups(ars(mxs,:)),'descend');
	    for ii = mxs-1:-1:1
		    [~,idx] = sort(arn(ii,ndx),'descend');
		    ndx = ndx(idx);
		    [~,idx] = sort(idn(ii,ndx),'descend');
		    ndx = ndx(idx);
		    [~,idx] = sort(nsGroups(ars(ii,ndx)),'descend');
		    ndx = ndx(idx);
	    end
    else
	    [~,ndx] = sort(ars(mxs,:)); % ascend
	    for ii = mxs-1:-1:1
		    [~,idx] = sort(arn(ii,ndx),'ascend');
		    ndx = ndx(idx);
		    [~,idx] = sort(idn(ii,ndx),'ascend');
		    ndx = ndx(idx);
		    [~,idx] = sort(ars(ii,ndx)); % ascend
		    ndx = ndx(idx);
	    end
    end
    %
    % Outputs %%
    %
    if ischar(A)
	    ndx = ndx(:);
	    B = A(ndx,:);
    else
	    ndx = reshape(ndx,size(A));
	    B = A(ndx);
    end
    %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsort
function grp = nsGroups(vec)
    % Groups in a cell array of char vectors, equivalent to [~,~,grp]=unique(vec);
    [vec,idx] = sort(vec);
    grp = cumsum([true(1,numel(vec)>0),~strcmp(vec(1:end-1),vec(2:end))]);
    grp(idx) = grp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsGroups
function nsChkRgx(rgx,chk)
    % Perform some basic sanity-checks on the supplied regular expression.
    chk = sprintf('^(%s)$',chk);
    assert(isempty(regexpi(rgx,chk,'once')),...
	    'SC:natsort:rgx:OptionMixUp',...
	    ['Second input <rgx> must be a regular expression that matches numbers.',...
	    '\nThe provided input "%s" looks like an optional argument (inputs 3+).'],rgx)
    if isempty(regexpi('0',rgx,'once'))
	    warning('SC:natsort:rgx:SanityCheck',...
		    ['Second input <rgx> must be a regular expression that matches numbers.',...
		    '\nThe provided regular expression does not match the digit "0", which\n',...
		    'may be acceptable (e.g. if literals, quantifiers, or lookarounds are used).'...
		    '\nThe provided regular expression: "%s"'],rgx)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsChkRgx
function nsAssert(txt,idx,eid,opt)
    % Throw an error if an option is overspecified.
    if nnz(idx)>1
	    error(sprintf('SC:natsort:%s:Overspecified',eid),...
		    ['The %s option may only be specified once.',...
		    '\nThe provided options:%s'],opt,sprintf(' "%s"',txt{idx}));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsAssert
function arr = ns1s2c(arr)
    % If scalar string then extract the character vector, otherwise data is unchanged.
    if isa(arr,'string') && isscalar(arr)
	    arr = arr{1};
    end
end
    
function [B,ndx,dbg] = natsortfiles(A,rgx,varargin)
    % Natural-order / alphanumeric sort of filenames or foldernames.
    %
    % (c) 2014-2023 Stephen Cobeldick
    %
    % Sorts text by character code and by number value. File/folder names, file
    % extensions, and path directories (if supplied) are sorted separately to
    % ensure that shorter names sort before longer names. For names without
    % file extensions (i.e. foldernames, or filenames without extensions) use
    % the 'noext' option. Use the 'xpath' option to ignore the filepath. Use
    % the 'rmdot' option to remove the folder names "." and ".." from the array.
    %
    %%% Example:
    % P = 'C:\SomeDir\SubDir';
    % S = dir(fullfile(P,'*.txt'));
    % S = natsortfiles(S);
    % for k = 1:numel(S)
    %     F = fullfile(P,S(k).name)
    % end
    %
    %%% Syntax:
    %  B = natsortfiles(A)
    %  B = natsortfiles(A,rgx)
    %  B = natsortfiles(A,rgx,<options>)
    % [B,ndx,dbg] = natsortfiles(A,...)
    %
    % To sort the elements of a string/cell array use NATSORT (File Exchange 34464)
    % To sort the rows of a string/cell/table use NATSORTROWS (File Exchange 47433)
    % To sort string/cells using custom sequences use ARBSORT (File Exchange 132263)
    %
    % File Dependency %%
    %
    % NATSORTFILES requires the function NATSORT (File Exchange 34464). Extra
    % optional arguments are passed directly to NATSORT. See NATSORT for case-
    % sensitivity, sort direction, number format matching, and other options.
    %
    % Explanation %%
    %
    % Using SORT on filenames will sort any of char(0:45), including the
    % printing characters ' !"#$%&''()*+,-', before the file extension
    % separator character '.'. Therefore NATSORTFILES splits the file-name
    % from the file-extension and sorts them separately. This ensures that
    % shorter names come before longer names (just like a dictionary):
    %
    % >> Af = {'test_new.m'; 'test-old.m'; 'test.m'};
    % >> sort(Af) % Note '-' sorts before '.':
    % ans =
    %     'test-old.m'
    %     'test.m'
    %     'test_new.m'
    % >> natsortfiles(Af) % Shorter names before longer:
    % ans =
    %     'test.m'
    %     'test-old.m'
    %     'test_new.m'
    %
    % Similarly the path separator character within filepaths can cause longer
    % directory names to sort before shorter ones, as char(0:46)<'/' and
    % char(0:91)<'\'. This example on a PC demonstrates why this matters:
    %
    % >> Ad = {'A1\B', 'A+/B', 'A/B1', 'A=/B', 'A\B0'};
    % >> sort(Ad)
    % ans =   'A+/B'  'A/B1'  'A1\B'  'A=/B'  'A\B0'
    % >> natsortfiles(Ad)
    % ans =   'A\B0'  'A/B1'  'A1\B'  'A+/B'  'A=/B'
    %
    % NATSORTFILES splits filepaths at each path separator character and sorts
    % every level of the directory hierarchy separately, ensuring that shorter
    % directory names sort before longer, regardless of the characters in the names.
    % On a PC separators are '/' and '\' characters, on Mac and Linux '/' only.
    %
    % Examples %%
    %
    % >> Aa = {'a2.txt', 'a10.txt', 'a1.txt'}
    % >> sort(Aa)
    % ans = 'a1.txt'  'a10.txt'  'a2.txt'
    % >> natsortfiles(Aa)
    % ans = 'a1.txt'  'a2.txt'  'a10.txt'
    %
    % >> Ab = {'test2.m'; 'test10-old.m'; 'test.m'; 'test10.m'; 'test1.m'};
    % >> sort(Ab) % Wrong number order:
    % ans =
    %    'test.m'
    %    'test1.m'
    %    'test10-old.m'
    %    'test10.m'
    %    'test2.m'
    % >> natsortfiles(Ab) % Shorter names before longer:
    % ans =
    %    'test.m'
    %    'test1.m'
    %    'test2.m'
    %    'test10.m'
    %    'test10-old.m'
    %
    %%% Directory Names:
    % >> Ac = {'A2-old\test.m';'A10\test.m';'A2\test.m';'A1\test.m';'A1-archive.zip'};
    % >> sort(Ac) % Wrong number order, and '-' sorts before '\':
    % ans =
    %    'A1-archive.zip'
    %    'A10\test.m'
    %    'A1\test.m'
    %    'A2-old\test.m'
    %    'A2\test.m'
    % >> natsortfiles(Ac) % Shorter names before longer:
    % ans =
    %    'A1\test.m'
    %    'A1-archive.zip'
    %    'A2\test.m'
    %    'A2-old\test.m'
    %    'A10\test.m'
    %
    % Input and Output Arguments %%
    %
    %%% Inputs (**=default):
    % A   = Array to be sorted. Can be the structure array returned by DIR,
    %       or a string array, or a cell array of character row vectors.
    % rgx = Optional regular expression to match number substrings.
    %     = []** uses the default regular expression (see NATSORT).
    % <options> can be supplied in any order:
    %     = 'rmdot' removes the dot directory names "." and ".." from the output.
    %     = 'noext' for foldernames, or filenames without filename extensions.
    %     = 'xpath' sorts by name only, excluding any preceding filepath.
    % Any remaining <options> are passed directly to NATSORT.
    %
    %%% Outputs:
    % B   = Array <A> sorted into natural sort order.      The same size as <A>.
    % ndx = NumericMatrix, generally such that B = A(ndx). The same size as <A>.
    % dbg = CellArray, each cell contains the debug cell array of one level
    %       of the filename/path parts, i.e. directory names, or filenames, or
    %       file extensions. Helps debug the regular expression (see NATSORT).
    %
    % See also SORT NATSORTFILES_TEST NATSORT NATSORTROWS ARBSORT IREGEXP
    % REGEXP DIR FILEPARTS FULLFILE NEXTNAME STRING CELLSTR SSCANF
    % Input Wrangling %%
    %
    fnh = @(c)cellfun('isclass',c,'char') & cellfun('size',c,1)<2 & cellfun('ndims',c)<3;
    %
    if isstruct(A)
	    assert(isfield(A,'name'),...
		    'SC:natsortfiles:A:StructMissingNameField',...
		    'If first input <A> is a struct then it must have field <name>.')
	    nmx = {A.name};
	    assert(all(fnh(nmx)),...
		    'SC:natsortfiles:A:NameFieldInvalidType',...
		    'First input <A> field <name> must contain only character row vectors.')
	    [fpt,fnm,fxt] = cellfun(@fileparts, nmx, 'UniformOutput',false);
	    if isfield(A,'folder')
		    fpt(:) = {A.folder};
		    assert(all(fnh(fpt)),...
			    'SC:natsortfiles:A:FolderFieldInvalidType',...
			    'First input <A> field <folder> must contain only character row vectors.')
	    end
    elseif iscell(A)
	    assert(all(fnh(A(:))),...
		    'SC:natsortfiles:A:CellContentInvalidType',...
		    'First input <A> cell array must contain only character row vectors.')
	    [fpt,fnm,fxt] = cellfun(@fileparts, A(:), 'UniformOutput',false);
	    nmx = strcat(fnm,fxt);
    elseif ischar(A)
	    assert(ndims(A)<3,...
		    'SC:natsortfiles:A:CharNotMatrix',...
		    'First input <A> if character class must be a matrix.') %#ok<ISMAT>
	    [fpt,fnm,fxt] = cellfun(@fileparts, num2cell(A,2), 'UniformOutput',false);
	    nmx = strcat(fnm,fxt);
    else
	    assert(isa(A,'string'),...
		    'SC:natsortfiles:A:InvalidType',...
		    'First input <A> must be a structure, a cell array, or a string array.');
	    [fpt,fnm,fxt] = cellfun(@fileparts, cellstr(A(:)), 'UniformOutput',false);
	    nmx = strcat(fnm,fxt);
    end
    %
    varargin = cellfun(@nsf1s2c, varargin, 'UniformOutput',false);
    ixv = fnh(varargin); % char
    txt = varargin(ixv); % char
    xtx = varargin(~ixv); % not
    %
    trd = strcmpi(txt,'rmdot');
    tnx = strcmpi(txt,'noext');
    txp = strcmpi(txt,'xpath');
    %
    nsfAssert(txt, trd, 'rmdot', '"." and ".." folder')
    nsfAssert(txt, tnx, 'noext', 'file-extension')
    nsfAssert(txt, txp, 'xpath', 'file-path')
    %
    chk = '(no|rm|x)(dot|ext|path)';
    %
    if nargin>1
	    nsfChkRgx(rgx,chk)
	    txt = [{rgx},txt(~(trd|tnx|txp))];
    end
    %
    % Path and Extension %%
    %
    % Path separator regular expression:
    if ispc()
	    psr = '[^/\\]+';
    else % Mac & Linux
	    psr = '[^/]+';
    end
    %
    if any(trd) % Remove "." and ".." dot directory names
	    ddx = strcmp(nmx,'.') | strcmp(nmx,'..');
	    fxt(ddx) = [];
	    fnm(ddx) = [];
	    fpt(ddx) = [];
	    nmx(ddx) = [];
    end
    %
    if any(tnx) % No file-extension
	    fnm = nmx;
	    fxt = [];
    end
    %
    if any(txp) % No file-path
	    mat = reshape(fnm,1,[]);
    else % Split path into {dir,subdir,subsubdir,...}:
	    spl = regexp(fpt(:),psr,'match');
	    nmn = 1+cellfun('length',spl(:));
	    mxn = max(nmn);
	    vec = 1:mxn;
	    mat = cell(mxn,numel(nmn));
	    mat(:) = {''};
	    %mat(mxn,:) = fnm(:); % old behavior
	    mat(permute(bsxfun(@eq,vec,nmn),[2,1])) =  fnm(:);  % TRANSPOSE bug loses type (R2013b)
	    mat(permute(bsxfun(@lt,vec,nmn),[2,1])) = [spl{:}]; % TRANSPOSE bug loses type (R2013b)
    end
    %
    if numel(fxt) % File-extension
	    mat(end+1,:) = fxt(:);
    end
    %
    % Sort Matrices %%
    %
    nmr = size(mat,1)*all(size(mat));
    dbg = cell(1,nmr);
    ndx = 1:numel(fnm);
    %
    for ii = nmr:-1:1
	    if nargout<3 % faster:
		    [~,idx] = natsort(mat(ii,ndx),txt{:},xtx{:});
	    else % for debugging:
		    [~,idx,gbd] = natsort(mat(ii,ndx),txt{:},xtx{:});
		    [~,idb] = sort(ndx);
		    dbg{ii} = gbd(idb,:);
	    end
	    ndx = ndx(idx);
    end
    %
    % Return the sorted input array and corresponding indices:
    %
    if any(trd)
	    tmp = find(~ddx);
	    ndx = tmp(ndx);
    end
    %
    ndx = ndx(:);
    %
    if ischar(A)
	    B = A(ndx,:);
    elseif any(trd)
	    xsz = size(A);
	    nsd = xsz~=1;
	    if nnz(nsd)==1 % vector
		    xsz(nsd) = numel(ndx);
		    ndx = reshape(ndx,xsz);
	    end
	    B = A(ndx);
    else
	    ndx = reshape(ndx,size(A));
	    B = A(ndx);
    end
    %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortfiles
function nsfChkRgx(rgx,chk)
    chk = sprintf('^(%s)$',chk);
    assert(~ischar(rgx)||isempty(regexpi(rgx,chk,'once')),...
	    'SC:natsortfiles:rgx:OptionMixUp',...
	    ['Second input <rgx> must be a regular expression that matches numbers.',...
	    '\nThe provided expression "%s" looks like an optional argument (inputs 3+).'],rgx)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsfChkRgx
function nsfAssert(txt,idx,eid,opt)
    % Throw an error if an option is overspecified.
    if nnz(idx)>1
	    error(sprintf('SC:natsortfiles:%s:Overspecified',eid),...
		    ['The %s option may only be specified once.',...
		    '\nThe provided options:%s'],opt,sprintf(' "%s"',txt{idx}));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nsfAssert
function arr = nsf1s2c(arr)
    % If scalar string then extract the character vector, otherwise data is unchanged.
    if isa(arr,'string') && isscalar(arr)
	    arr = arr{1};
    end
end

% From “Fast and Robust Curve Intersections” v2.0.0.0 on MathWorks File Exchange

function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
    %INTERSECTIONS Intersections of curves.
    %   Computes the (x,y) locations where two curves intersect.  The curves
    %   can be broken with NaNs or have vertical segments.
    %
    % Example:
    %   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
    %
    % where X1 and Y1 are equal-length vectors of at least two points and
    % represent curve 1.  Similarly, X2 and Y2 represent curve 2.
    % X0 and Y0 are column vectors containing the points at which the two
    % curves intersect.
    %
    % ROBUST (optional) set to 1 or true means to use a slight variation of the
    % algorithm that might return duplicates of some intersection points, and
    % then remove those duplicates.  The default is true, but since the
    % algorithm is slightly slower you can set it to false if you know that
    % your curves don't intersect at any segment boundaries.  Also, the robust
    % version properly handles parallel and overlapping segments.
    %
    % The algorithm can return two additional vectors that indicate which
    % segment pairs contain intersections and where they are:
    %
    %   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
    %
    % For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
    % (how far along this segment the intersection is).  For example, if I(k) =
    % 45.25 then the intersection lies a quarter of the way between the line
    % segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
    % the vector J and the segments in (X2,Y2).
    %
    % You can also get intersections of a curve with itself.  Simply pass in
    % only one curve, i.e.,
    %
    %   [X0,Y0] = intersections(X1,Y1,ROBUST);
    %
    % where, as before, ROBUST is optional.
    
    % Version: 2.0, 25 May 2017
    % Author:  Douglas M. Schwarz
    % Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
    % Real_email = regexprep(Email,{'=','*'},{'@','.'})
    
    
    % Theory of operation:
    %
    % Given two line segments, L1 and L2,
    %
    %   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
    %   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
    %
    % we can write four equations with four unknowns and then solve them.  The
    % four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
    % L1 and L2, t1 is the distance from the starting point of L1 to the
    % intersection relative to the length of L1 and t2 is the distance from the
    % starting point of L2 to the intersection relative to the length of L2.
    %
    % So, the four equations are
    %
    %    (x1(2) - x1(1))*t1 = x0 - x1(1)
    %    (x2(2) - x2(1))*t2 = x0 - x2(1)
    %    (y1(2) - y1(1))*t1 = y0 - y1(1)
    %    (y2(2) - y2(1))*t2 = y0 - y2(1)
    %
    % Rearranging and writing in matrix form,
    %
    %  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
    %        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
    %   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
    %        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
    %
    % Let's call that A*T = B.  We can solve for T with T = A\B.
    %
    % Once we have our solution we just have to look at t1 and t2 to determine
    % whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
    % line segments cross and we can include (x0,y0) in the output.
    %
    % In principle, we have to perform this computation on every pair of line
    % segments in the input data.  This can be quite a large number of pairs so
    % we will reduce it by doing a simple preliminary check to eliminate line
    % segment pairs that could not possibly cross.  The check is to look at the
    % smallest enclosing rectangles (with sides parallel to the axes) for each
    % line segment pair and see if they overlap.  If they do then we have to
    % compute t1 and t2 (via the A\B computation) to see if the line segments
    % cross, but if they don't then the line segments cannot cross.  In a
    % typical application, this technique will eliminate most of the potential
    % line segment pairs.
    
    
    % Input checks.
    if verLessThan('matlab','7.13')
	    error(nargchk(2,5,nargin)) %#ok<NCHKN>
    else
	    narginchk(2,5)
    end
    
    % Adjustments based on number of arguments.
    switch nargin
	    case 2
		    robust = true;
		    x2 = x1;
		    y2 = y1;
		    self_intersect = true;
	    case 3
		    robust = x2;
		    x2 = x1;
		    y2 = y1;
		    self_intersect = true;
	    case 4
		    robust = true;
		    self_intersect = false;
	    case 5
		    self_intersect = false;
    end
    
    % x1 and y1 must be vectors with same number of points (at least 2).
    if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
		    length(x1) ~= length(y1)
	    error('X1 and Y1 must be equal-length vectors of at least 2 points.')
    end
    % x2 and y2 must be vectors with same number of points (at least 2).
    if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
		    length(x2) ~= length(y2)
	    error('X2 and Y2 must be equal-length vectors of at least 2 points.')
    end
    
    
    % Force all inputs to be column vectors.
    x1 = x1(:);
    y1 = y1(:);
    x2 = x2(:);
    y2 = y2(:);
    
    % Compute number of line segments in each curve and some differences we'll
    % need later.
    n1 = length(x1) - 1;
    n2 = length(x2) - 1;
    xy1 = [x1 y1];
    xy2 = [x2 y2];
    dxy1 = diff(xy1);
    dxy2 = diff(xy2);
    
    
    % Determine the combinations of i and j where the rectangle enclosing the
    % i'th line segment of curve 1 overlaps with the rectangle enclosing the
    % j'th line segment of curve 2.
    
    % Original method that works in old MATLAB versions, but is slower than
    % using binary singleton expansion (explicit or implicit).
    % [i,j] = find( ...
    % 	repmat(mvmin(x1),1,n2) <= repmat(mvmax(x2).',n1,1) & ...
    % 	repmat(mvmax(x1),1,n2) >= repmat(mvmin(x2).',n1,1) & ...
    % 	repmat(mvmin(y1),1,n2) <= repmat(mvmax(y2).',n1,1) & ...
    % 	repmat(mvmax(y1),1,n2) >= repmat(mvmin(y2).',n1,1));
    
    % Select an algorithm based on MATLAB version and number of line
    % segments in each curve.  We want to avoid forming large matrices for
    % large numbers of line segments.  If the matrices are not too large,
    % choose the best method available for the MATLAB version.
    if n1 > 1000 || n2 > 1000 || verLessThan('matlab','7.4')
	    % Determine which curve has the most line segments.
	    if n1 >= n2
		    % Curve 1 has more segments, loop over segments of curve 2.
		    ijc = cell(1,n2);
		    min_x1 = mvmin(x1);
		    max_x1 = mvmax(x1);
		    min_y1 = mvmin(y1);
		    max_y1 = mvmax(y1);
		    for k = 1:n2
			    k1 = k + 1;
			    ijc{k} = find( ...
				    min_x1 <= max(x2(k),x2(k1)) & max_x1 >= min(x2(k),x2(k1)) & ...
				    min_y1 <= max(y2(k),y2(k1)) & max_y1 >= min(y2(k),y2(k1)));
			    ijc{k}(:,2) = k;
		    end
		    ij = vertcat(ijc{:});
		    i = ij(:,1);
		    j = ij(:,2);
	    else
		    % Curve 2 has more segments, loop over segments of curve 1.
		    ijc = cell(1,n1);
		    min_x2 = mvmin(x2);
		    max_x2 = mvmax(x2);
		    min_y2 = mvmin(y2);
		    max_y2 = mvmax(y2);
		    for k = 1:n1
			    k1 = k + 1;
			    ijc{k}(:,2) = find( ...
				    min_x2 <= max(x1(k),x1(k1)) & max_x2 >= min(x1(k),x1(k1)) & ...
				    min_y2 <= max(y1(k),y1(k1)) & max_y2 >= min(y1(k),y1(k1)));
			    ijc{k}(:,1) = k;
		    end
		    ij = vertcat(ijc{:});
		    i = ij(:,1);
		    j = ij(:,2);
	    end
	    
    elseif verLessThan('matlab','9.1')
	    % Use bsxfun.
	    [i,j] = find( ...
		    bsxfun(@le,mvmin(x1),mvmax(x2).') & ...
		    bsxfun(@ge,mvmax(x1),mvmin(x2).') & ...
		    bsxfun(@le,mvmin(y1),mvmax(y2).') & ...
		    bsxfun(@ge,mvmax(y1),mvmin(y2).'));
	    
    else
	    % Use implicit expansion.
	    [i,j] = find( ...
		    mvmin(x1) <= mvmax(x2).' & mvmax(x1) >= mvmin(x2).' & ...
		    mvmin(y1) <= mvmax(y2).' & mvmax(y1) >= mvmin(y2).');
	    
    end
    
    
    % Find segments pairs which have at least one vertex = NaN and remove them.
    % This line is a fast way of finding such segment pairs.  We take
    % advantage of the fact that NaNs propagate through calculations, in
    % particular subtraction (in the calculation of dxy1 and dxy2, which we
    % need anyway) and addition.
    % At the same time we can remove redundant combinations of i and j in the
    % case of finding intersections of a line with itself.
    if self_intersect
	    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
    else
	    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
    end
    i(remove) = [];
    j(remove) = [];
    
    % Initialize matrices.  We'll put the T's and B's in matrices and use them
    % one column at a time.  AA is a 3-D extension of A where we'll use one
    % plane at a time.
    n = length(i);
    T = zeros(4,n);
    AA = zeros(4,4,n);
    AA([1 2],3,:) = -1;
    AA([3 4],4,:) = -1;
    AA([1 3],1,:) = dxy1(i,:).';
    AA([2 4],2,:) = dxy2(j,:).';
    B = -[x1(i) x2(j) y1(i) y2(j)].';
    
    % Loop through possibilities.  Trap singularity warning and then use
    % lastwarn to see if that plane of AA is near singular.  Process any such
    % segment pairs to determine if they are colinear (overlap) or merely
    % parallel.  That test consists of checking to see if one of the endpoints
    % of the curve 2 segment lies on the curve 1 segment.  This is done by
    % checking the cross product
    %
    %   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
    %
    % If this is close to zero then the segments overlap.
    
    % If the robust option is false then we assume no two segment pairs are
    % parallel and just go ahead and do the computation.  If A is ever singular
    % a warning will appear.  This is faster and obviously you should use it
    % only when you know you will never have overlapping or parallel segment
    % pairs.
    
    if robust
	    overlap = false(n,1);
	    warning_state = warning('off','MATLAB:singularMatrix');
	    % Use try-catch to guarantee original warning state is restored.
	    try
		    lastwarn('')
		    for k = 1:n
			    T(:,k) = AA(:,:,k)\B(:,k);
			    [unused,last_warn] = lastwarn; %#ok<ASGLU>
			    lastwarn('')
			    if strcmp(last_warn,'MATLAB:singularMatrix')
				    % Force in_range(k) to be false.
				    T(1,k) = NaN;
				    % Determine if these segments overlap or are just parallel.
				    overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			    end
		    end
		    warning(warning_state)
	    catch err
		    warning(warning_state)
		    rethrow(err)
	    end
	    % Find where t1 and t2 are between 0 and 1 and return the corresponding
	    % x0 and y0 values.
	    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
	    % For overlapping segment pairs the algorithm will return an
	    % intersection point that is at the center of the overlapping region.
	    if any(overlap)
		    ia = i(overlap);
		    ja = j(overlap);
		    % set x0 and y0 to middle of overlapping region.
		    T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			    min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		    T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			    min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		    selected = in_range | overlap;
	    else
		    selected = in_range;
	    end
	    xy0 = T(3:4,selected).';
	    
	    % Remove duplicate intersection points.
	    [xy0,index] = unique(xy0,'rows');
	    x0 = xy0(:,1);
	    y0 = xy0(:,2);
	    
	    % Compute how far along each line segment the intersections are.
	    if nargout > 2
		    sel_index = find(selected);
		    sel = sel_index(index);
		    iout = i(sel) + T(1,sel).';
		    jout = j(sel) + T(2,sel).';
	    end
    else % non-robust option
	    for k = 1:n
		    [L,U] = lu(AA(:,:,k));
		    T(:,k) = U\(L\B(:,k));
	    end
	    
	    % Find where t1 and t2 are between 0 and 1 and return the corresponding
	    % x0 and y0 values.
	    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
	    x0 = T(3,in_range).';
	    y0 = T(4,in_range).';
	    
	    % Compute how far along each line segment the intersections are.
	    if nargout > 2
		    iout = i(in_range) + T(1,in_range).';
		    jout = j(in_range) + T(2,in_range).';
	    end
    end
end

% Plot the results (useful for debugging).
% plot(x1,y1,x2,y2,x0,y0,'ok');

function y = mvmin(x)
    % Faster implementation of movmin(x,k) when k = 1.
    y = min(x(1:end-1),x(2:end));
end

function y = mvmax(x)
    % Faster implementation of movmax(x,k) when k = 1.
    y = max(x(1:end-1),x(2:end));
end


