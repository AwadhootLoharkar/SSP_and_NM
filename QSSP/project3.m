%Excercise 4.8.5
mystartdefaults;

%units
recipunit = 1.0e+10; %angstrom
ekinscale =((hbar*recipunit)^2/(2*elm))/qel; %kinetic energy convertion

%input data
a = 4.08; %Lattice spacing [Angstrom]: value for Au
nband = 4; %number of bands to be plotted

%reciprocal lattice vectors
ngx = 10; %minimal number of rows/columns of Hamiltonian which corresponds to the number of reciprocal lattice vectors
ngx = floor(ngx/2); % rounds ngx/2 to the integer value below it
fprintf('%d x %d Hamiltonian\n',2*ngx+1,2*ngx+1); %Hamiltonian should be NxN for N an odd number

ng = 2*ngx+1; %we need the double values since we will define one for G1 but also one for -G1
n=0;
%defining the G for plus and minus values and the zero value (eq. 4.20)
for i=-ng:ng
    n=n+1;
    G(n)=i; %reciprocal lattice vectors, we are working in 2pi/a units (integers) to define G
end

%Defining the Fourier components of potential (V_G)
V_G = zeros(1,2*ng+1);
n=0;
for i=-ng:ng
    n=n+1;
    if (abs(i)==1)
        V_G(n)=1; % 1 eV gap
    end
    if (abs(i)==2)
        V_G(n)=1.5; % 1.5 eV gap
    end
end

%sampling recirpocal space (discretizing it) in 2pi/a units
k_max = 1.5; %maximum value in the horizontal axis
k = -k_max:0.01:k_max;
nk = length(k); %number of k points

%Drawing BZ boundaries
E_max_plot = 35;
cb = E_max_plot *(-1).^(abs(round(k)));

%Hamiltonian (for flat potential Vg = 0)
H = zeros(ngx,ngx);%initialization of the 11x11 Hamiltonian

%Hamiltonian (for Vg != 0) - non-diagonal elements of the Hamiltonian
% In this case we set the non-diagonal elements before the loop, since
% these elements do not depend on k and won't be affected in the loop, the
% loop just affects the diagonal (kinetic energy) elements
j = 0;
for jx=-ngx:ngx
    j=j+1;
    i=0;
    for ix=-ngx:ngx
        i=i+1;
        H(i,j)=V_G(ng+1+jx-ix);
    end
end

%Main loop
for m=1:nk
    %Kinetic energy (diagonal elements)
    n=0;
    for i=-ngx:ngx
        n=n+1;
        H(n,n)=ekinscale*(2*pi/a)^2*(k(m)-G(ng+1+i))^2; %equation 4.21; we multiply by (2*pi/a)^2 to recover the units 
    end

    %Diagonalization
    if(any(abs(H-H'))>1e-12)
        error('Hamiltonian is not Hermitian');
    end
    [v,ev] = eig(H); %storing the vectors v and eigenvectors ev
    E = real(diag(ev)); %making sure that the eigenvalues are real numbers
    [E,perm] = sort(E); %ordering eigenvalues from smaller to greather
    %v=v(:,perm); %for Vg =! 0 later

    %storing data as colums of array bands (1:nk, 1:nband)
    for i=1:nband %(recall that nband=4)
        bands(m,i)=E(i); %energy bands given by the eigenvalues
    end
end


plot(k,bands);hold on;
xlim([-1.5 1.5]);
ylim([0 35]);
plot(k,cb,'color','black');

