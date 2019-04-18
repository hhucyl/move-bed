clear
clc
filenum = [1:99];
prefix_name = {'test_periodic1_'};
N = ceil(2.0*pi*20);
alpha = 2*pi/N;
Rh = 20;
L = [];
color = {'k*','r*'};
for i = 1:numel(filenum)
    name = strcat(char(prefix_name),num2str(filenum(i),'%04d'),'.h5');
    Fh = h5read(char(name),char('/PForceh'));
    V = h5read(char(name),char('/PVeloc'));
%     plot(i,Fh(1),'k*')
    FFh(i,1) = Fh(1);
%     hold on
    drawnow
    Nx = h5read(char(name),char('/Nx'));
    Box = [0 Nx-1];
    pos = h5read(char(name),char('/Pposition'));
    ppos(i,1)=pos(1);
    if(pos(1)<Box(1)+Rh || pos(1)>Box(2)-Rh)
        plot(i,abs(Fh(1)),char(color(2)));
        hold on
    else
        plot(i,abs(Fh(1)),char(color(1)));
        hold on
    end
    
end
ylabel('\itF_y')
xlabel('Step')
fontsize = 20;
set(gca,'FontName','Times New Roman', 'FontSize',fontsize)
set(get(gca,'XLabel'),'FontName','Times New Roman', 'FontSize',fontsize)
set(get(gca,'YLabel'),'FontName','Times New Roman', 'FontSize',fontsize)