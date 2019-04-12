clear
clc
prefix_name = {'/home/pzhang/chen/move-bed/'};
file_name = {'test_move1_0020.h5'};
kkx = 200:600;
kky = 500:1600;

name = strcat(prefix_name,file_name);
vel = h5read(char(name),char('/Velocity_0'));
nx  = h5read(char(name),char('/Nx'));
ny  = h5read(char(name),char('/Ny'));
U   = vel(1:3:end-2);
V   = vel(3:3:end);
u   = reshape(U,[nx,ny]);
v   = reshape(V,[nx,ny]);
kkx = 1:nx;
kky = 1:ny;
uu  = u(kkx,kky)';
vv  = v(kkx,kky)';
[M,N] = size(uu);
[x,y] = meshgrid(1:N,1:M);
cav = curl(x,y,uu,vv);
pcolor(x,y,cav)
shading interp
colorbar
% caxis('manual')
% caxis(gca,[1e-3 -7e-3])
caxis([-0.003 0.002])
hold on
pos = h5read(char(name),'/Pposition');
viscircles([pos(1),pos(2)],10);
axis equal

