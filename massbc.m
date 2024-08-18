function [pl,ql,pr,qr] = massbc(~,~,~,~,~)
%Boundary conditions, no diffusion at the walls
    pl=0;
    ql=1;               
    pr=0;
    qr=1;               
end