clear
clc
prefix_name = {'/home/pzhang/chen/move-bed/'};
prefix_name = strcat(char(prefix_name),'test_mvbedrw_1_');
num = [0:999];
R = 10*0.8;
layer = [189,149,109,69];
color = {'k','b','g','y'};
for i = 1:numel(num)
    
    filename = strcat(prefix_name,num2str(num(i),'%04d'),'.h5');
    pos = h5read(filename,'/RWPposition');
    ppos = h5read(filename,'/Pposition');
    vel = h5read(filename,'/Velocity_0');
    ga = h5read(filename,'/Gamma');
    nx = h5read(filename,'/Nx');
    ny = h5read(filename,'/Ny');
    Velx = reshape(vel(1:3:end-2),[nx,ny]);
    Vely = reshape(vel(2:3:end-1),[nx,ny]);
    Ga = reshape(ga,[nx,ny]);
    Nrwp = numel(pos)/3;
    Pos = [pos(1:3:end-2),pos(2:3:end-1),pos(3:3:end)];
    Np = numel(ppos)/6;
    Ppos = [ppos(1:3:(3*Np-2)),ppos(2:3:(3*Np-1))];
    RR = zeros(1,Np)+R;
    subplot(2,3,[1,2,3])
    viscircles(Ppos,RR);
    hold on
    plot(Pos(:,1),Pos(:,2),'b.')
    for j=1:numel(layer)
        plot(0:nx-1,ones(nx).*layer(j),char(strcat(color(j),'--')),'linewidth',2)
    end
    axis equal
    axis([0 nx-1,0 ny-1])
    subplot(2,3,[4])
    for j=1:numel(layer)
        plot(0:nx-1,Vely(:,layer(j)),char(strcat(color(j),'-')),'linewidth',1)
        hold on
    end
    subplot(2,3,[5])
    for j=1:numel(layer)
        yy = Vely(:,layer(j));
        plot(0:nx-1,(yy-min(yy))./(max(yy)-min(yy)),char(strcat(color(j),'-')),'linewidth',1)
        hold on
    end
    subplot(2,3,[6])
    vvely = Vely(:,layer);
    for j=1:numel(layer)
        temp = vvely(:,j);
        flux(i,j) = sum(temp(temp>0));
        plot(flux(:,j),char(strcat(color(j),'*')))
        hold on
    end
    
    drawnow
    jpgname = strcat(prefix_name,num2str(num(i),'%04d'),'.jpg');
    saveas(gcf,jpgname)
    clf
end