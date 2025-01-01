%% PHYSICAL CONSTANTS MOST FREQUENTLY USED 
%% IN PHOTONICS, QUANTUM SOLID STATE PHYSICS & NANOTECHNOLOGIES

%  Link to complete list of updated SI fundamental physical constants:
%  https://physics.nist.gov/cuu/Constants/Table/allascii.txt

%% THE SEVEN EXACT DEFINING CONSTANTS OF THE SI UNIT SYSTEM (2019 UPDATE)

hyperfine = 9192631770;      % Hyperfine transition frequency of Cs_133 [Hz]
celeritas = 299792458;       % Speed of light in vacuum [m/s]
Planck    = 6.62607015E-34;  % Planck constant [Js] 
qel       = 1.602176634E-19; % Elementary charge [C]
kB        = 1.380640E-23;    % Boltzmann constant [J/K]
Avogadro  = 6.02214076E23;   % Avogadro constant [1/mole]
kcd       = 683;             % Luminous efficacy [candela=lumen/Watt]
          % in detecting 540 THz radiation (Green light @ 555.016 nm).
          % Maximum possible luminous efficacy by a photopic observer.
          % Originally, peak sensitivity of "average" human eye.

%% PHYSICAL CONSTANTS OF ELECTROMAGNETISM

epsilon0  = 8.8541878128E-12;% Electrical constant [F/m]
                             % (formerly vacuum dielectric permittivity) 
mu0       = 1.25663706212E-6;% Magnetic constant [N/A^2]
                             % (formerly vacuum magnetic permeability)
                             % close to 4*pi*1.E-7 =1.25663706143E-6
Klitzing  = 25812.80745;     % von Klitzing's constant = Planck/qel^2 [Ohm]
                             % respecting significant digits 

%% UNITS USED IN QUANTUM PHYSICS

hbar      = 1.054571817E-34;  % Reduced Planck's constant = Planck/(2*pi)
                              % respecting significant digits [Js]
Angstroem = 1.0E-10;          % Angström [m]
Dalton    = 1.66053906660E-27;% Atomic mass unit [kg] = 1 Dalton
                              % = mass of Carbon_12 atom / 12
elm       = 9.1093837015E-31; % Electron mass [kg]
nem       = 1.67492749804E-27;% Neutron mass [kg]
prm       = 1.67262192369E-27;% Proton mass [kg]
elecint   = qel^2/(4*pi*epsilon0);% Electron-electron interaction scale [N*m^2]
Bohr      = hbar^2/(elm*elecint)/Angstroem;     % Bohr radius [Angström]
Rydberg   = (elm/(2*hbar^2)) * (elecint)^2/qel; % Rydberg [eV]
Hartree   = 2*Rydberg;                          % Hartree [eV]