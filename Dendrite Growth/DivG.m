% Computing divergence of gradient

function y=DivG(g,f,dx,dy)
cs=@circshift;

g_1p=cs(g,-1,1);
g_1n=cs(g,1,1);

g_2p=cs(g,-1,2);
g_2n=cs(g,1,2);

f_1p=cs(f,-1,1);
f_1n=cs(f,1,1);
         
f_2p=cs(f,-1,2);
f_2n=cs(f,1,2);
        
y1=0.5*((f_1p+f).*(g_1p-g)-(f_1n+f).*(g-g_1n))/dx^2;
y2=0.5*((f_2p+f).*(g_2p-g)-(f_2n+f).*(g-g_2n))/dy^2;
y=y1+y2;
y(1,:)=3*y(2,:)-3*y(3,:)+y(4,:);
y(end,:)=3*y(end-1,:)-3*y(end-2,:)+y(end-3,:);
yb=y([1,end],[1,end]);
y(:,1)=3*y(:,2)-3*y(:,3)+y(:,4);
y(:,end)=3*y(:,end-1)-3*y(:,end-2)+y(:,end-3);
y([1,end],[1,end])=0.5*(yb+y([1,end],[1,end]));
end
