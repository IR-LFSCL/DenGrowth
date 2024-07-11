% Dendrite Growth
% Zn: component 1; Zn ion: component 2; Na ion: component 3;
% Ac ion: component 4

%% Simulation parameters
time0=clock();

% Basic constants
D0=3.0e-10;                                      % m^2/s
w0=1.0e7;                                        % J/m^3
eps0=4.0e-7;                                     % J/m
F0=96485.3383;                                   % C/mol
R0=8.314;                                        % J/(K*mol)
T0=298;                                          % K

% Characteristic quantities
d_c=sqrt(eps0/w0);                               % length
t_c=eps0/(w0*D0);                                % time
f_c=w0;                                          % energy density
C_c=w0/(R0*T0);                                  % molar concentration
pot_c=R0*T0/F0;                                  % electric potential
sig_c=w0*D0*(F0/R0/T0)^2;                        % conductivity
j_c=D0*sqrt(w0/eps0)*w0/pot_c;                   % current density
v_c=d_c/t_c;                                     % velocity

% Physical parameters
De=3.68e-15;                                     % m^2/s
Ds=D0;                                           % m^2/s
se=1.0e7;                                        % S/m
ss=1.0;                                          % S/m
Ce=1.09e5;                                       % mol/m^3
Cs=5.56e4;                                       % mol/m^3
v_int=1.0e-6;                                    % m/s
pot0=-0.5;                                       % V
j0=0.2;                                          % C/(m^2*s)

alpha=0.5;
delta=0.02;
omega=0.05;
m0=6;
z=[0,2,1,-1,0];

% Grid and Output parameters
Nx=200;
Ny=200;
dx=1.0;
dy=1.0;
dt=0.001;
Nt=50000;
Sp=10000;

err=1.0e-6;

% Reduced quantities
De_r=De/D0;
Ds_r=Ds/D0;
se_r=se/sig_c;
ss_r=ss/sig_c;

Ce_r=Ce/C_c;
Cs_r=Cs/C_c;
dV=1/Ce_r-1/Cs_r;

pote_r=-1.0/pot_c;
pots_r=0.0/pot_c;
pot0_r=pot0/pot_c;
j_r=j0/j_c;

v_r=v_int/v_c;

field={'x1','x2','x3','x4','pot'};
field_p='x1';
folder_d='Data';

% x mol/L ZnAc2, y mol/L NaAc
C2=@(x,y)1000*x;
C3=@(x,y)1000*y;
C4=@(x,y)2000*x+1000*y;
C5=@(x,y)Cs-3000*x-2000*y;

x_min=0.001;
x2s=C2(0.5,0.5)/Cs*(1-x_min);
x3s=C3(0.5,0.5)/Cs*(1-x_min);
x4s=C4(0.5,0.5)/Cs*(1-x_min);
x5s=C5(0.5,0.5)/Cs*(1-x_min);

%% Defining needed field variables
x1=zeros(Nx,Ny);
x2=zeros(Nx,Ny);
x3=zeros(Nx,Ny);
x4=zeros(Nx,Ny);
x5=zeros(Nx,Ny);
pot=zeros(Nx,Ny);

%% Assigning initial condition
for i=1:Nx
    for j=1:Ny
        if i<=15
            x1(i,j)=1-4*x_min;
            x2(i,j)=x_min;
            x3(i,j)=x_min;
            x4(i,j)=x_min;
            x5(i,j)=x_min;
            pot(i,j)=pote_r;
        else
            x1(i,j)=x_min;
            x2(i,j)=x2s;
            x3(i,j)=x3s;
            x4(i,j)=x4s;
            x5(i,j)=x5s;
            pot(i,j)=pots_r;
        end
    end
end

DensityPlot(field_p,[0,1],0);
OutputData(folder_d,field,0);
       
%% Updating field variables and outputing data
for n=1:Nt
    % Computing derivatives of phi concerning x and y
    Dx_x1=Dx(x1,dx);
    Dy_x1=Dy(x1,dy);
    Dphi=sqrt((Dx_x1*dx).^2+(Dy_x1*dy).^2);
    
    % Computing auxiliary quantities
    theta=atan2(Dy_x1,Dx_x1);
    cos_v=cos(m0*theta);
    eps=(1+delta*cos_v).^2;
    Deps=-2*delta*m0*(1+delta*cos_v).*sin(m0*theta);
    
    Mm=(x1*De_r+(1-x1)*Ds_r);
    sm=x1*se_r+(1-x1)*ss_r;
    Cm=1./(x1/Ce_r+(1-x1)/Cs_r);
    Vmol=x1/Ce_r+(1-x1)/Cs_r;
    Gmol=x1.*log(x1)+x2.*log(x2)+x3.*log(x3)+x4.*log(x4)+...
        x5.*log(x5)+(z(1)*x1+z(2)*x2+z(3)*x3+z(4)*x4+...
        z(5)*x5).*pot;
    
    % Updating pot and assighing its boundary condition
    pot=EllSol(pot,sm,0,dx,dy,err);
    
    Dx1_f=Cm.*(log(x1)+1+z(1)*pot)-dV*Cm.^2.*Gmol+...
        x1.*(1-x1).*(1-2*x1)-DivG(phi,eps,dx,dy)-...
        Dy(0.5*Deps.*Dx_phi,dy)+Dx(0.5*Deps.*Dy_phi,dx);
    Dx2_f=Cm.*(log(x2)+1+z(2)*pot);
    Dx3_f=Cm.*(log(x3)+1+z(3)*pot);
    Dx4_f=Cm.*(log(x4)+1+z(4)*pot);
    Dx5_f=Cm.*(log(x5)+1+z(5)*pot);
    
    DivJ1=DivG(Dx2_f-Dx1_f,Mm*x1.*x2,dx,dy)+...
        DivG(Dx3_f-Dx1_f,Mm*x1.*x3,dx,dy)+...
        DivG(Dx4_f-Dx1_f,Mm*x1.*x4,dx,dy)+...
        DivG(Dx5_f-Dx1_f,Mm*x1.*x5,dx,dy);
    DivJ2=DivG(Dx1_f-Dx2_f,Mm*x2.*x1,dx,dy)+...
        DivG(Dx3_f-Dx2_f,Mm*x2.*x3,dx,dy)+...
        DivG(Dx4_f-Dx2_f,Mm*x2.*x4,dx,dy)+...
        DivG(Dx5_f-Dx2_f,Mm*x2.*x5,dx,dy);
    DivJ3=DivG(Dx1_f-Dx3_f,Mm*x3.*x1,dx,dy)+...
        DivG(Dx2_f-Dx3_f,Mm*x3.*x2,dx,dy)+...
        DivG(Dx4_f-Dx3_f,Mm*x3.*x4,dx,dy)+...
        DivG(Dx5_f-Dx3_f,Mm*x3.*x5,dx,dy);
    DivJ4=DivG(Dx1_f-Dx4_f,Mm*x4.*x1,dx,dy)+...
        DivG(Dx2_f-Dx4_f,Mm*x4.*x2,dx,dy)+...
        DivG(Dx3_f-Dx4_f,Mm*x4.*x3,dx,dy)+...
        DivG(Dx5_f-Dx4_f,Mm*x4.*x5,dx,dy);
    
    eta=(Dx2_f-Dx1_f)./(z(2)*Cm);
    jRac=j_r*(exp(-alpha*z(2)*eta)-exp((1-alpha)*z(2)*eta));
    Rate=Dphi.*jRac;
    
    
    % Computing derivatives of phi and con concerning t
    Ran=omega*2*(rand(Nx,Ny)-0.5).*x1.*(1-x1);
    Dt_x1=Vmol.*(-DivJ1+Rate)+Ran;
    Dt_x2=Vmol.*(-DivJ2-Rate);
    Dt_x3=Vmol.*(-DivJ3);
    Dt_x4=Vmol.*(-DivJ4);
    
    % Updating phi and assighing its boundary condition
    x1=x1+dt*Dt_x1;
    x1(x1<0)=0;
    x1(x1>1)=1;
    x2=x2+dt*Dt_x2;
    x2(x2<0)=0;
    x2(x2>1)=2;
    x3=x3+dt*Dt_x3;
    x3(x3<0)=0;
    x3(x3>1)=1;
    x4=x4+dt*Dt_x4;
    x4(x4<0)=0;
    x4(x4>1)=1;
    x5=1-x1-x2-x3-x4;
    x5(x5<0)=0;
    x5(x5>1)=1;
    
    x_sum=x1+x2+x3+x4+x5;
    x1=x1./x_sum;
    x2=x2./x_sum;
    x3=x3./x_sum;
    x4=x4./x_sum;
    x5=1-x1-x2-x3-x4;
    
    % Updating con and assighing its boundary condition
%     con=con+dt*Dt_con;
%     con(con<0.0)=0.0;
%     con=BC_con(con,phi);
    
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
