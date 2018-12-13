clear
clc
prefix = '/home/pzhang/chen/move-bed/';
middle = 'test_mvbed_';
for i = 1:409
    name = strcat(prefix,middle,num2str(i,'%04d'),'.h5');
    nx = h5read(name,'/Nx');
    ny = h5read(name,'/Ny');
    vel = h5read(name,'/Velocity_0');
    p = h5read(name,'/Density_0');
    ga = h5read(name,'/Gamma');
    ii = 1:nx*ny;
    vx = reshape(vel(3*(ii-1)+1),[nx,ny]);
    vy = reshape(vel(3*(ii-1)+1),[nx,ny]);
    P = reshape(p,[nx,ny]);
    Ga = reshape(ga,[nx,ny]);
    subplot(221)
    title(i)
    quiver(vx,vy)
    subplot(222)
    vvx = sum(vx,1);
    plot(-vvx,1:ny)
    subplot(223);
    plot(P(:,401))
    Pa(:,i) = P(:,401);
    xlim([1,400])
    subplot(224)
    vvx1 = vx(200,:);
    pv = P(200,:).*vvx1.^2;
    
    Y(i,1) = sum(pv);
    plot(i,Y(i,1),'b*')
    hold on
    drawnow
end
figure
plot(-vvx./numel(1:ny),1:ny)
xlabel('\itV')
ylabel('depth')
figure()
sp = size(Pa);
Paa = sum(Pa,2);
plot((Paa-min(Paa))./(max(Paa)-min(Paa)))
ylabel('\itp^*')
xlim([0 400])