%% DEFINE THE BAND STRUCTURE PLOT PARAMETERS
semiconductor = 'ZnTe'; % Choose the semiconductor
bs_step = 7; % Choose the step for the BS scattering plot

%% READ THE BAND STRUCTURE FROM A FILE
materials = {'Si', 'Ge', 'Sn', 'GaP', 'GaAs', 'AlSb', 'InP', 'GaSb', ...
    'InAs', 'InSb', 'ZnS', 'ZnSe', 'ZnTe', 'CdTe', 'Empty lattice'};
m = find(strcmp(materials, semiconductor));
filename = strcat(int2str(m), 'bandstructure.dat');
fid = fopen(filename, 'r');
nqpath = fscanf(fid, '%d', 1); % Read nqpath from the file
nband = fscanf(fid, '%d', 1); % Read nband from the file
Eband = zeros(nband, nqpath); % Initialize Eband with the correct size
for ii = 1:nqpath
    for jj = 1:nband
        qpath(5, ii) = fscanf(fid, '%f', 1);
        Eband(jj, ii) = fscanf(fid, '%f', 1);
    end
    fscanf(fid, '%c', 1);
end
tix = fscanf(fid, '%f', [1, inf]); % Read tix
til = {}; % Initialize til as an empty cell array
line = fgetl(fid); % Read the first line
index = 1; % Initialize index
while ischar(line) % Continue until fgetl returns -1
    if ~isempty(line) % If the line is not empty
        if mod(index, 2) == 0 % If index is even
            line = sprintf('\n\n%s', line); % Add two new lines before the line
        end
        til{end + 1} = line; % Add the line to til
        index = index + 1; % Increment index
    end
    line = fgetl(fid); % Read the next line
end
fclose(fid);

%% PLOT THE BAND STRUCTURE
title_pos = [[1.3, 5]; [1.3, 5]; [1.3, 5]; [1.3, 5]; [1.3, 5]; [1.3, 5]; [1.3, 5]; ...
    [1.3, 5]; [1.3, 5]; [1.3, 5]; [1.3, 8]; [1.3, 8]; [1.4, 7]; [1.4, 6]];
ylimits = [[-6, 7]; [-5, 7]; [-4, 6]; [-4, 7]; [-4, 7]; [-4, 7]; [-4, 7]; ...
    [-3, 6.4]; [-4, 7]; [-3, 6]; [-3, 10]; [-3, 9]; [-3, 9]; [-4, 8]];
figure;
hold on;
plot(qpath(5, :), Eband, '-', 'color', 'black', 'linewidth', 1.3); % Plot the band structure with a line
plot(qpath(5, 1:bs_step:end), Eband(:, 1:bs_step:end), 'o', 'color', ...
    'black', 'MarkerSize', 5, 'linewidth', 1.3, 'MarkerFaceColor', 'white'); % Plot the band structure with circles and step 'bs_step'
hold off

% Define custom tick labels and positions
custom_labels = {'L','[-1,-1,-1]', '\Gamma','[1,0,0]', 'X','[-1,-1,0]', 'K','[-1,-1,0]', '\Gamma'};
custom_positions = tix; % Define the custom x-axis positions
title(materials{m}, 'fontsize', 26, 'position', title_pos(m, :));
ylabel('E   (eV)', 'FontSize', 24);
xlim([0, qpath(5, nqpath)]);
ylim(ylimits(m, :));
set(gca, 'xtick', custom_positions);
set(gca, 'xticklabel', custom_labels, 'FontSize', 24);
set(gca, 'ytick', ylimits(m, 1):1:ylimits(m, 2));
set(gca, 'TickLength', [0.03, 0.02]);
set(gca, 'Box', 'on');
set(gca, 'linewidth', 1);
