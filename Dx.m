% Computing first derivative concerning x

function y=Dx(g,dx)
cs=@circshift;

g_1p=cs(g,-1,1);
g_1n=cs(g,1,1);
    
y=(g_1p-g_1n)/(2*dx);
y(1,:)=3*y(2,:)-3*y(3,:)+y(4,:);
y(end,:)=3*y(end-1,:)-3*y(end-2,:)+y(end-3,:);
end

