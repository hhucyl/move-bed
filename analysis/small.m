clear
clc
prefix = {'/media/pzhang/My Book/dune_shape/small/1/'};
fname = {'test_mvbed_c_'};
num = 1:1019;
index = 800;
for i = 1:numel(num)
    name = strcat(prefix,fname,num2str(num(i),'%04d'),'.h5');
    nx = h5read(char(name),char('/Nx'));
    ny = h5read(char(name),char('/Ny'));
    v = h5read(char(name),char('/Velocity_0'));
    vx = reshape(v(1:3:end-2),[nx,ny]);
    vy = reshape(v(2:3:end-1),[nx,ny]);
    V(i,1) = sum(sum(sqrt(vx(:,index:end).^2+vy(:,index:end).^2)))/double(numel(index:ny)*nx);
    if(i>931)
        plot(10*(i-931)+931,V(i,1),'b*')
    else
        plot(i,V(i,1),'b*')
    end
    hold on
    drawnow
end