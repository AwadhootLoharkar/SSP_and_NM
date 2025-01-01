%% INPUT PARAMETERS FOR COMPUTING
%% THE ELECTRON BAND STRUCTURE OF FCC SEMICONDUCTORS
%% USING EMPIRICAL PSEUDOPOTENTIALS 

semiconductor='Si'

%% ADJUSTING ZERO ENERGY LEVEL
%  Accepted strings (case insensitive):
%  'E_F'   or 'Fermi level'
%  'E_vac' or 'vacuum level'
%  'bcb'   or 'bottom of conduction band'
%  'tvb'   or 'top of valence band' 

adjust_zero='top of valence band' 

%% Tolerance for recognizing zero

tol=1e-12; 

%% PARAMETERS IN RECIPROCAL SPACE
cutoff=21   % Deal only with |G|^2 < cutoff [2*pi/spacing]^2
            % in Hamiltonian

Gs_max=11   % |G|^2 of highest non zero Fourier coefficients in
            % expanding potential [2*pi/spacing]^2

%% FOR COMPUTING THE DISPERSION RELATIONS
%% SAMPLING OF BRILLOUIN ZONE = PATH ASSEMBLING SEGMENTS

if dispersion_relation
    nbands=16    % Number of bands to be stored in output file
    BZstep=0.02 % Step along path in BZ
    
    qs(1:3,1)= [0.5 0.5 0.5]';    qs_str{1} = 'L';      % start point 1
    qe(1:3,1)= [0 0 0]';          qe_str{1} = '\Gamma'; % end   point 1
    ticksep(1)=0.25; 
    
    qs(1:3,2)= [0 0 0]';          qs_str{2} = '\Gamma'; % start point 2
    qe(1:3,2)= [1 0 0]';          qe_str{2} = 'X';      % end   point 2
    ticksep(2)=0.2; 
    
    qs(1:3,3)= [1 1 0]';          qs_str{3} = 'X';      % start point 3
    qe(1:3,3)= [0.75 0.75 0]';    qe_str{3} = 'K';      % end   point 3
    ticksep(3)=0.5;
    
    qs(1:3,4)= [0.75 0.75 0]';    qs_str{4} = 'K';      % start point 4
    qe(1:3,4)= [0 0 0]';          qe_str{4} = '\Gamma'; % end   point 4
    ticksep(4)=0.25;
end

%% DATA OF EMPIRICAL PSEUDOPOTENTIALS
%% FOR COMPUTING ELECTRON ENERGY BANDS
%% OF 14 FCC SEMICONDUCTORS
%
%  from M.L. Cohen & T.K. Bergstresser, 
%       Phys. Rev. vol.141, p.789 (1966)

materials={'Si';'Ge';'Sn';'GaP';'GaAs';'AlSb';'InP';'GaSb'; ...
           'InAs';'InSb';'ZnS';'ZnSe';'ZnTe';'CdTe';'Empty lattice'};

%% MATERIAL IDENTIFIER

m = find(strcmp(materials,semiconductor));

if isempty(m)
    error('Material not recognized.')
end

%% DIRECT LATTICE UNIT VECTORS AND ATOMIC POSITIONS
%% OF FCC SEMICONDUCTORS

a(1:3,1)= [0.5 0.5 0.0]' ;    % Direct lattice unit vector 1
a(1:3,2)= [0.0 0.5 0.5]' ;    % Direct lattice unit vector 2
a(1:3,3)= [0.5 0.0 0.5]' ;    % direct lattice unit vector 3
cell_volume = a(:,1)' * cross(a(:,2),a(:,3)); %  = det(a)

tau(1:3,1) = [ 0.125  0.125  0.125]' ; % Position of atom 1 in primitive cell
tau(1:3,2) = [-0.125 -0.125 -0.125]' ; % Position of atom 2 in primitive cell

%% LATTICE SPACINGS [Angström]

ls(1:14)= [5.43 5.66 6.49 5.44 5.64 6.13 5.86 6.12 6.04 6.48 5.41 ...
           5.65 6.07 6.41];

%% PSEUDOPOTENTIAL FORM FACTORS [Rydberg]

% Remark: ff(i,1) = V0 = V_{G=0} for material i
%                 = constant adjusting zero of potential
%                   could be set in this table

%%            V0  VS3  VS8  VS11 VA3  VA4  VA11
ff(1,:) = [ 0.00 -0.21 0.04 0.08 0.00 0.00 0.00]; % Si
ff(2,:) = [ 0.00 -0.23 0.01 0.06 0.00 0.00 0.00]; % Ge
ff(3,:) = [ 0.00 -0.20 0.00 0.04 0.00 0.00 0.00]; % Sn
ff(4,:) = [ 0.00 -0.22 0.03 0.07 0.12 0.07 0.02]; % GaP
ff(5,:) = [ 0.00 -0.23 0.01 0.06 0.07 0.05 0.01]; % GaAs
ff(6,:) = [ 0.00 -0.21 0.02 0.06 0.06 0.04 0.02]; % AlSb
ff(7,:) = [ 0.00 -0.23 0.01 0.06 0.07 0.05 0.01]; % InP
ff(8,:) = [ 0.00 -0.22 0.00 0.05 0.06 0.05 0.01]; % GaSb
ff(9,:) = [ 0.00 -0.22 0.00 0.05 0.08 0.05 0.03]; % InAs
ff(10,:)= [ 0.00 -0.20 0.00 0.04 0.06 0.05 0.01]; % InSb
ff(11,:)= [ 0.00 -0.22 0.03 0.07 0.24 0.14 0.04]; % ZnS
ff(12,:)= [ 0.00 -0.23 0.01 0.06 0.18 0.12 0.03]; % ZnSe
ff(13,:)= [ 0.00 -0.22 0.00 0.05 0.13 0.10 0.01]; % ZnTe
ff(14,:)= [ 0.00 -0.20 0.00 0.04 0.15 0.09 0.04]; % CdTe

%% WORK FUNCTIONS [eV]
% Data from various references.
% Uncertainty may be large (0.1 to 1 eV).
% Sources of uncertainty: doping;
%                         quality of crystal;
%                         quality of sample surface.

wf(1)  = 4.85; % Si
wf(2)  = 4.75; % Ge
wf(3)  = 4.42; % Sn
wf(4)  = 4.34; % GaP
wf(5)  = 4.69; % GaAS
wf(6)  = 4.46; % AlSb
wf(7)  = 4.65; % InP
wf(8)  = 4.45; % GaSb  (4.1 to 4.8 eV)
wf(9)  = 4.95; % InAs
wf(10) = 4.57; % InSb
wf(11) = 7.00; % ZnS  (7 to 8 eV)
wf(12) = 5.11; % ZnSe
wf(13) = 5.80; % ZnTe (5.3 to 5.8 eV)
wf(14) = 5.70; % CdTe
