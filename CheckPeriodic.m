clear
clc
filenum = [1:99];
prefix_name = {'test_periodic1_'};
N = ceil(2.0*pi*20);
alpha = 2*pi/N;
Rh = 20;

for i = 1:numel(filenum)
    name = strcat(char(prefix_name),num2str(filenum(i),'%04d'),'.h5');
    Fh = h5read(char(name),char('/PForceh'));
    V = h5read(char(name),char('/PVeloc'))
    plot(i,Fh(2),'k*')
    FFh(i,1) = Fh(1);
    hold on
    drawnow

    pos = h5read(char(name),char('/Pposition'));
    ppos(i,1)=pos(1);
    ppos(i,2)=pos(4);
end