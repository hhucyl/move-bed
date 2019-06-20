clear
clc
prefix = '/home/pzhang/chen/move-bed/';
prefix = '/media/pzhang/My Book/move-bed-tmp/move_bed_2/';
prefix = '/media/pzhang/My Book/move-bed-tmp/macondo/0.001r_20.0Ga_0.3gap/';
for iiii = 1:1
% if(iiii==1)
%     middle = 'test_mvbed_';
% else
%     middle = strcat('test_mvbed_',num2str(iiii-1),'_');
    middle = 'test_mvbed_c_';
% end
R = 10;
LP = 40;
num = 1:986;
for i = 1:numel(num)
    name = strcat(prefix,middle,num2str(num(i),'%04d'),'.h5');
    nx = h5read(name,'/Nx');
    ny = h5read(name,'/Ny');
    vel = h5read(name,'/Velocity_0');
    p = h5read(name,'/Density_0');
    ga = h5read(name,'/Gamma');
    ppos = h5read(name,'/Pposition');
    pvec = h5read(name,'/PVeloc');
    NP = numel(ppos)/6;
    ip = 1:NP;
    Ppos = [ppos(3*(ip-1)+1),ppos(3*(ip-1)+2),ppos(3*(ip-1)+3)];
    Lxp = Ppos(401:440,1);
    Lyp = Ppos(401:440,2);
    Lyp_n = (Lyp-min(Lyp))./(max(Lyp)-min(Lyp));
    [llxp(:,i),ki] = sort(Lxp);
    llyp(:,i) = Lyp(ki);
    ii = 1:nx*ny;
    vx = reshape(vel(3*(ii-1)+1),[nx,ny]);
    vy = reshape(vel(3*(ii-1)+1),[nx,ny]);
    P = reshape(p,[nx,ny]);
    Ga = reshape(ga,[nx,ny]);
    gas = sum(Ga,1);
    kga = find(gas<1);
%     subplot(221)
%     title(i)
%     quiver(vx,vy)
    subplot(231)
    vvx = mean(vx,1);
    plot(-vvx,1:ny)
    xlabel('V_a')
    ylabel('depth')
    subplot(232);
    ppp = P(:,kga(2));
    ppp_n = (ppp-min(ppp))./(max(ppp)-min(ppp));
    plot(ppp_n)
    hold on 
    plot(Lxp,Lyp_n,'rs')
    Pa(:,i) = ppp_n;
    ylabel('p^*')
    xlim([1,nx])
    hold off
    subplot(233)
    plot(1:nx,mean(P,2))
    subplot(234)
    vvx1 = vx(0.5*nx,:);
    pv = P(0.5*nx,:).*vvx1;
    Y(i,1) = -sum(pv);
    plot(i+(iiii-1)*numel(num),Y(i,1),'b*')
    ylabel('\rhov^2')
    hold on
    subplot(235)
    dp = sum(P(1,:))-sum(P(nx,:));
    plot(i+(iiii-1)*numel(num),dp,'b*')
    hold on
    ylabel('\DeltaP')
%     subplot(335)
%     [temp,kpp] = sort(Ppos(:,2),'descend');
%     for j=1:LP
% %         viscircles(Ppos(400+j,1:2),R,'color','r');
%         plot(Ppos(400+j,1),Ppos(400+j,2),'s')
%         hold on
%     end
%     hold off
% %     axis equal
    subplot(236)
    ppv(i,1) = (sum(pvec(1:3*NP).^2))^0.5;
    plot(i+(iiii-1)*numel(num),ppv(i,1),'b*')
    hold on
    ylabel('\v_pa^2')
    drawnow
end

% figure()
% sp = size(Pa);
% Paa = sum(Pa,2);
% plot((Paa-min(Paa))./(max(Paa)-min(Paa)))
% ylabel('\itp^*')
% lax = mean(llxp,2);
% lay = mean(llyp,2);
% lay_n = (lay-min(lay))./(max(lay)-min(lay));
% hold on
% plot(1:20:nx-40,lay_n,'rs')
end