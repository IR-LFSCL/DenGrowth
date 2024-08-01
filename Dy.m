% Computing first derivative concerning y

function y=Dy(g,dy)
cs=@circshift;

g_2p=cs(g,-1,2);
g_2n=cs(g,1,2);
    
y=(g_2p-g_2n)/(2*dy);
y(:,1)=3*y(:,2)-3*y(:,3)+y(:,4);
y(:,end)=3*y(:,end-1)-3*y(:,end-2)+y(:,end-3);
end
