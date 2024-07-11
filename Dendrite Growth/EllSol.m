% Elliptic Solver

function y=EllSol(g,f,s,dx,dy,err)
cs=@circshift;

Niter=10000;
ryx=dy^2/dx^2;

phi=evalin('base','phi');

f_1p=cs(f,-1,1);
f_1n=cs(f,1,1);
         
f_2p=cs(f,-1,2);
f_2n=cs(f,1,2);

f1=0.5*(f_1p+f_1n+2*f);
f2=0.5*(f_2p+f_2n+2*f)*ryx;

for iter=1:Niter
    g_old=g;
    
    g_1p=cs(g,-1,1);
    g_1n=cs(g,1,1);
    
    g_2p=cs(g,-1,2);
    g_2n=cs(g,1,2);
    
    gf1=0.5*((f_1p+f).*g_1p+(f_1n+f).*g_1n);
    gf2=0.5*((f_2p+f).*g_2p+(f_2n+f).*g_2n)*ryx;
    
    g=(gf1+gf2-s*dx^2)./(f1+f2);
    g=BC_pot(g,phi);
    
    err_max=max(abs(g-g_old),[],'all');
    if err_max<=err || iter==Niter
        y=g;
        break;
    end
end
end
