% Dendrite Growth

%% Simulation parameters
time0=clock();

L0=2.5e-6;             % m^3/(J*s)
w0=1.3344e7;           % J/m^3
eps0=4.17e-7;          % J/m
F0=96485.3383;         % C/mol

Lc=0.00206;            % 1/s
RT=8.314*300;          % J/mol
De=3.68e-15;           % m^2/s
Ds=3.68e-10;           % m^2/s
se=6.3e7;              % S/m
ss=2.67;               % S/m
Ce=9.71e4;             % mol/m^3
Cs=1.1e4;              % mol/m^3
V0=1.0;                % V

d0=sqrt(eps0/w0);      % m
t0=1/(L0*w0);          % s

Nx=200;
Ny=200;
dx=0.8;
dy=0.8;
dt=0.0001;
Nt=50000;
Sp=10000;

alpha=0.5;
delta=0.02;
omega=7.5;
m0=6;
z0=1;
err=1.0e-6;

Lc_r=Lc*t0;
F0_r=z0*F0/RT;

Ds_r=Ds*t0/d0^2;
De_r=De*t0/d0^2;

s0=z0*F0*Ce*L0*eps0/V0;
se_r=se/s0;
ss_r=ss/s0;
ce_r=0.0;
cs_r=1.0;
ve_r=-0.7;
vs_r=0.0;
C_sca=Ce/Cs;
k_pot=(vs_r-ve_r)/(Nx-1);

field={'phi','con','pot','con_Na','con_Zn'};
field_p='phi';
folder_d='Data';

h=@(x)x.^3.*(6*x.^2-15*x+10);
Dh=@(x)30*x.^2.*(1-x).^2;

%% Defining needed field variables
phi=zeros(Nx,Ny);
con=zeros(Nx,Ny);
pot=zeros(Nx,Ny);
con_Na=zeros(Nx,Ny);
con_Zn=zeros(Nx,Ny);

%% Assigning initial condition
for i=1:Nx
    for j=1:Ny
        if i<=15
            phi(i,j)=1.0;
            pot(i,j)=-0.7;
        else
            con(i,j)=1.0;
        end
    end
end

DensityPlot(field_p,[0,1],0);
OutputData(folder_d,field,0);
       
%% Updating field variables and outputing data
for n=1:Nt
    % Computing derivatives of phi concerning x and y
    Dx_phi=Dx(phi,dx);
    Dy_phi=Dy(phi,dy);
    
    % Computing auxiliary quantities
    theta=atan2(Dy_phi,Dx_phi);
    cos_v=cos(m0*theta);
    eps=(1+delta*cos_v).^2;
    Deps=-2*delta*m0*(1+delta*cos_v).*sin(m0*theta);
    
    he=h(phi);
    Dm=he*De_r+(1-he)*Ds_r;
    sm=he*se_r+(1-he)*ss_r;
    
    % Updating pot and assighing its boundary condition
    Dhe=Dh(phi);
    pot=EllSol(pot,sm,0,dx,dy,err);
    
    % Computing derivatives of phi and con concerning t
    Dt_phi=Dy(0.5*Deps.*Dx_phi,dy)-Dx(0.5*Deps.*Dy_phi,dx)+...
        DivG(phi,eps,dx,dy)-phi.*(1-phi).*(1-2*phi)+...
        omega*2*(rand(Nx,Ny)-0.5).*Dhe-Lc_r*Dhe.*...
        (exp((1-alpha)*F0_r*pot)-con.*exp(-alpha*F0_r*pot));
    Dt_con=DivG(con,Dm,dx,dy)+DivG(pot,F0_r*Dm.*con,dx,dy)-...
        C_sca*Dt_phi;
    
    % Updating phi and assighing its boundary condition
    phi=phi+dt*Dt_phi;
    phi=BC_phi(phi);
    
    % Updating con and assighing its boundary condition
    con=con+dt*Dt_con;
    con(con<0.0)=0.0;
    con(con>1.0)=1.0;
    con=BC_con(con,phi);
    
    con_Na=0.1.*con;
    con_Zn=0.9.*con;
    
    % Outputing data at specified time level
    if mod(n,Sp)==0
        DensityPlot(field_p,[0,1],n);
        OutputData(folder_d,field,n);
        fprintf('done steps: %d\n',n);
    end
end

%% Displaying consumed time
con_time = etime(clock(),time0);
fprintf('Time Consumed: %10d\n',con_time);
