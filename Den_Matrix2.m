% Dendrite Growth

%% Simulation parameters
time0=clock();

L0=1.5e-6;             % m^3/(J*s)
w0=1.3344e7;           % J/m^3
eps0=4.17e-7;          % J/m
F0=96485.3383;         % C/mol

Lc=0.00206;            % 1/s
RT=8.314*300;          % J/mol
De=3.68e-13;           % m^2/s
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
Niter=10000;
ryx=dy^2/dx^2;

alpha=0.5;
delta=0.02;
omega=0.05;
m0=6;
zZn=2;
zNa=1;
err=1.0e-6;
pot0=-0.34;

Lc_r=Lc*t0;
F0_r=zZn*F0/RT;

% Ds_r=Ds*t0/d0^2;
% De_r=De*t0/d0^2;
Ds_r=840;
De_r=0.84;

s0=zZn*F0*Ce*L0*eps0/V0;
se_r=se/s0;
ss_r=ss/s0;

C_sca=Ce/Cs;
CZn0=0.7;
CNa0=1-CZn0;

field={'phi','cZn','cNa','pot'};
field_p='phi';
folder_d='Data';

h=@(x)x.^3.*(6*x.^2-15*x+10);
Dh=@(x)30*x.^2.*(1-x).^2;
cs=@circshift;

%% Defining needed field variables
phi=zeros(Nx,Ny);
cZn=zeros(Nx,Ny);
cNa=zeros(Nx,Ny);
pot=zeros(Nx,Ny);

%% Assigning initial condition
for i=1:Nx
    for j=1:Ny
        if i<=10
            phi(i,j)=1.0;
            pot(i,j)=pot0;
        end
        cZn(i,j)=CZn0*(1-phi(i,j));
        cNa(i,j)=CNa0*(1-phi(i,j));
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
    Dt_phi1=Dy(0.5*Deps.*Dx_phi,dy)-Dx(0.5*Deps.*Dy_phi,dx)+...
        DivG(phi,eps,dx,dy)-phi.*(1-phi).*(1-2*phi)+...
        omega*2*(rand(Nx,Ny)-0.5).*Dhe;
    
    sm_1p=cs(sm,-1,1);
    sm_1n=cs(sm,1,1); 
    sm_2p=cs(sm,-1,2);
    sm_2n=cs(sm,1,2);
    
    we1=0.5*(sm_1p+sm_1n+2*sm);
    we2=0.5*(sm_2p+sm_2n+2*sm)*ryx;
    
    for iter=1:Niter
        pot_old=pot;
        
        pot_1p=cs(pot,-1,1);
        pot_1n=cs(pot,1,1);
        pot_2p=cs(pot,-1,2);
        pot_2n=cs(pot,1,2);
    
        pot_we1=0.5*((sm_1p+sm).*pot_1p+(sm_1n+sm).*pot_1n);
        pot_we2=0.5*((sm_2p+sm).*pot_2p+(sm_2n+sm).*pot_2n)*ryx;
        
        Dt_phi=Dt_phi1-Lc_r*Dhe.*(exp((1-alpha)*F0_r*...
            pot)-(cZn/CZn0).*exp(-alpha*F0_r*pot));
        pot=(pot_we1+pot_we2-Dt_phi*dx^2)./(we1+we2);
        pot(1,:)=pot0;
        pot(end,:)=0.0;
        pot(:,1)=pot(:,2);
        pot(:,end)=pot(:,end-1);
    
        err_max=max(abs(pot-pot_old),[],'all');
        if err_max<=err || iter==Niter
            break;
        end
    end
    
    % Computing derivatives of phi and con concerning t
    Dt_phi=Dt_phi1-Lc_r*Dhe.*(exp((1-alpha)*F0_r*...
        pot)-(cZn/CZn0).*exp(-alpha*F0_r*pot));
    Dt_cZn=DivG(cZn,Dm,dx,dy)+DivG(pot,F0_r*Dm.*cZn,dx,dy)-...
        C_sca*Dt_phi;
    Dt_cNa=DivG(cNa,Dm,dx,dy)+DivG(pot,F0_r*Dm.*cNa,dx,dy);
    
    % Updating phi and assighing its boundary condition
    phi=phi+dt*Dt_phi;
    phi(1,:)=1;
    phi(end,:)=0;
    phi(:,1)=phi(:,2);
    phi(:,end)=phi(:,end-1);
    
    % Updating cZn/cNa and assighing its boundary condition 
    cZn=cZn+dt*Dt_cZn;
    cZn(cZn<0.0)=0.0;
    cZn(cZn>CZn0)=CZn0;
    cZn(1,:)=0.0;
    cZn(end,:)=CZn0;
    cZn(:,1)=cZn(:,2);
    cZn(:,end)=cZn(:,end-1);
    
    cNa=cNa.*(1-phi);
    cNa=cNa+dt*Dt_cNa;
    cNa(1,:)=0.0;
    cNa(end,:)=CNa0;
    cNa(:,1)=cNa(:,2);
    cNa(:,end)=cNa(:,end-1);
    
    % Outputing data at specified time level
    if mod(n,Sp)==0
        DensityPlot(field_p,[0,1],n);
        OutputData(folder_d,field,n);
        
        mainf=cd;
        ss=num2str(n);
        fig_name=['phi','_',ss];
        surf(phi,phi,'EdgeColor','none')
        axis equal off; caxis([0,1]); colorbar; view([0,90]);
        set(gcf,'PaperPositionMode','manual',...
            'PaperUnits','points','PaperPosition',...
            [0 0 360 360],'Name',fig_name);
        cd([mainf,'/','Atlas']);
        print(gcf,'-dtiff',fig_name);
        cd(mainf);
        
        fprintf('done steps: %d\n',n);
    end
end

%% Displaying consumed time
con_time = etime(clock(),time0);
fprintf('Time Consumed: %10d\n',con_time);
