%% constants, lenth in nm, energy in eV
N=500;          % # of grid points

E_f=0.01;       % fermi energy in eV
E_g=0.1;        % band gap in eV
V_ds=0.1;       % drain-source voltage in V
V_g=0;          % gate potential Psi_g=-e*V_g in eV
d_ox=1;         % oxide thickness in nm
d_ch=10;        % channel thickness in nm
e=util.const.e; % elementary charge 
k_0=8.85e-12;   % dielectric constant
k_Si=11.68;     % dielectric constant for Silicon
k_ox=3.9;       % dielectric constant oxide
geo=1;          % Geometriefaktor # of gates, 'w' wrapgate
l_ch=10;        % channel length
l_ds=10;        % length of drain and source regions
N_dot=0;        % # Dopands


lambda=sqrt(k_Si/k_ox*d_ch*d_ox/geo); 

%% create sparse laplacian
tic;
super=zeros(N,1);
super(2:end)=1;
super(2)=2;
sub=zeros(N,1);
sub(1:end-1)=1;
sub(end-1)=2;
middle=zeros(N,1);
middle(:)=-2;
L_sparse=spdiags([super,middle,sub],[1,0,-1],N,N);
toc;

%% calculate source, drain and channel regions
% calc overall length
l_g=2*l_ds+l_ch;
% calculate spacing
a=l_g/(N-1);
% calculate indices
n_left=floor(l_ds/a);   % source region 1 to n-left
n_right=N-n_left;       % drain region n_right-end

%% assign Psi_g and Psi_b for all regions
%create empty arrays
Psi_g=zeros(N,1);
Psi_bi=zeros(N,1);

% soure both 0

% channel
Psi_g(n_left+1:n_right-1) =-V_g;
Psi_bi(n_left+1:n_right-1)=E_f+E_g/2.;

%drain
Psi_g(n_right:end) =0;
Psi_bi(n_right:end)=-V_ds;

% create charge density
rho=zeros(N,1);

%% solve eq. for Psi_f
% solve sparse (much quicker!)
tic;
Psi_f=(L_sparse-spdiags(zeros(N,1)+1/lambda^2,0,N,N))\((rho+N_dot)/k_0/k_Si-1/lambda^2*(Psi_g+Psi_bi));
toc;
%% Plot
figure, plot((0:N-1).*a,Psi_f);


