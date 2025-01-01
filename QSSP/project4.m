%% Adding paths for necessary functions and setting up defaults
addpath("../../start-kit-student");
addpath("../bin");

mystartdefaults; % Load default parameters and constants
tolerance = 1E-12; % Numerical tolerance for hermitian check

%% Constants and Conversion Factors
reciprocalUnit = 1E10; % Reciprocal length unit [m^-1], scale of the Brillouin zone
% Conversion factor: (ħ * k)^2 / (2 * m_e) -> energy in electron volts (eV)
energyScaleFactor = (hbar * reciprocalUnit)^2 / (2 * elm * qel);

%% Input Parameters for the Model
latticeSpacing = 5; % Lattice spacing [Å] (angstroms)
numReciprocalVectors = 10; % Number of reciprocal lattice vectors to consider
numKPoints = 100; % Number of k-points to sample in the Brillouin zone
kPointMax = 1.5; % Maximum value of k on the x-axis in units of [2π/a]
numBands = 6; % Number of energy bands to compute and plot

%% Ensure the number of reciprocal vectors is odd for matrix symmetry
numReciprocalVectors = floor(numReciprocalVectors / 2) * 2 + 1;
% Generate G vectors (reciprocal lattice vectors) centered at 0
G_vectors = -floor(numReciprocalVectors / 2):floor(numReciprocalVectors / 2);

% Check if the number of bands to store is valid
if (numBands > numReciprocalVectors)
    error("The number of bands must be less than or equal to the number of G vectors.");
end

%% Initialize k-points and eigenvalue storage
% k-points range from -kPointMax to +kPointMax
kPoints = linspace(-kPointMax, kPointMax, numKPoints);
eigenValues = zeros(numBands, numKPoints); % Matrix to store the computed eigenvalues

%% Loop over each k-point to compute the band structure
for kIndex = 1:length(kPoints)
    currentK = kPoints(kIndex); % Current k-point in the loop

    %% Initialize Hamiltonian Matrix
    Hamiltonian = zeros(numReciprocalVectors, numReciprocalVectors); % Hamiltonian matrix for this k-point

    %% Construct the Diagonal Elements of the Hamiltonian (Kinetic Energy)
    % Using the form E = ħ^2 * (k - G)^2 / (2m)
    for i = 1:length(G_vectors)
        gVector = G_vectors(i);
        kineticEnergy = energyScaleFactor * (currentK - gVector)^2;
        Hamiltonian(i, i) = kineticEnergy; % Assign to diagonal element
    end

    %% Check if Hamiltonian is Hermitian (H == H^†)
    if (any(abs(Hamiltonian - Hamiltonian') > tolerance))
        error("Hamiltonian is not Hermitian. Check matrix construction.");
    end

    %% Diagonalize the Hamiltonian to find eigenvalues (energies) and eigenvectors
    [eigenVectors, eigenMatrix] = eig(Hamiltonian); 
    eigenValuesList = real(diag(eigenMatrix)); % Extract eigenvalues (energies)

    %% Sort eigenvalues and store the lowest numBands bands
    [sortedEigenValues, sortOrder] = sort(eigenValuesList);
    eigenValues(:, kIndex) = sortedEigenValues(1:numBands); % Store first numBands eigenvalues
end

%% Plotting the Electron Band Structure
figure('Name', 'Electron Band Structure - Flat Potential', 'NumberTitle', 'off');
hold on;
% Plot each energy band
for bandIndex = 1:numBands
    plot(kPoints, eigenValues(bandIndex, :), 'DisplayName', ['Band ', num2str(bandIndex)]);
end
xlabel("Wavevector k [2\pi/a]");
ylabel("Energy E(k) [eV]");
xlim([-kPointMax, kPointMax]); % Set x-axis limits
ylim([0, max(max(eigenValues))]); % Set y-axis limits based on max energy
title('Electron Band Structure for a 1D Crystal with Flat Potential');
legend('show'); % Show legend for different bands
grid on;
movegui('northwest'); % Position the figure window
