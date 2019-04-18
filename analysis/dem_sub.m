clear
clc
prefix = '/home/pzhang/chen/move-bed/';
middle = 'test_move';
color = {'b*','k*','r*'};
for iii=1:3
for i = 1:84
    name = strcat(prefix,middle,num2str(iii),'_',num2str(i,'%04d'),'.h5');
    pos = h5read(name,'/Pposition');
    nx = h5read(name,'/Nx');
    ny = h5read(name,'/Ny');
    Np = numel(pos)/6;
    ii = 1:Np;
    pv = h5read(name,'/PVeloc');
    PV = [pv(3*(ii-1)+1),pv(3*(ii-1)+2),pv(3*(ii-1)+3)];
    plot(i,PV(2),char(color(iii)))
    pvv(i,iii) = PV(2);
    hold on
    drawnow
end
end