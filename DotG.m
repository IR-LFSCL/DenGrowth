% Computing dot product of gradients

function y=DotG(g,f,dx,dy)
cs=@circshift;

g_1p=cs(g,-1,1);
g_1n=cs(g,1,1);

g_2p=cs(g,-1,2);
g_2n=cs(g,1,2);

f_1p=cs(f,-1,1);
f_1n=cs(f,1,1);
         
f_2p=cs(f,-1,2);
f_2n=cs(f,1,2);
        
y1=0.25*(f_1p-f_1n).*(g_1p-g_1n)/dx^2;
y2=0.25*(f_2p-f_2n).*(g_2p-g_2n)/dy^2;
y=y1+y2;
y(1,:)=3*y(2,:)-3*y(3,:)+y(4,:);
y(end,:)=3*y(end-1,:)-3*y(end-2,:)+y(end-3,:);
yb=y([1,end],[1,end]);
y(:,1)=3*y(:,2)-3*y(:,3)+y(:,4);
y(:,end)=3*y(:,end-1)-3*y(:,end-2)+y(:,end-3);
y([1,end],[1,end])=0.5*(yb+y([1,end],[1,end]));
end
