clear
clc
N = ceil(2.0*pi*20);
alpha = 2*pi/N;
Rh = 20;
X = [75,25];
L = 70;
GX = [5,25];
nx = 70;
ny = 50;
% for i=0:nx-1
%     for j=0:ny-1
%         plot(i,j,'gs')
%         hold on
%         drawnow
%     end
% end
depth = 0;
for i=0:N-1;
    r = [Rh*cos(i*alpha)+X(1),Rh*sin(i*alpha)+X(2)];
    gr=[Rh*cos(i*alpha)+GX(1),Rh*sin(i*alpha)+GX(2)];
    ixs = max(floor(r(1) - depth),0.0);
    ixe = min(ceil(r(1) + depth),nx-1);
    iys = max(floor(r(2) - depth),0.0);
    iye = min(ceil(r(2) + depth),ny-1);
    plot(r(1),r(2),'b*')
    hold on
    plot([r(1),X(1)],[r(2),X(2)],'-b')
    xx = [ixs ixe;iys iye];
    for ix = ixs:ixe
        for iy = iys:iye
            plot(ix,iy,'bo')
        end
    end
    
    ixsg = max(floor(gr(1) - depth),0.0);
    ixeg = min(ceil(gr(1) + depth),nx-1);
    iysg = max(floor(gr(2) - depth),0.0);
    iyeg = min(ceil(gr(2) + depth),ny-1);
    gxx = [ixsg ixeg;iysg iyeg];
    plot(gr(1),gr(2),'k*')
    plot([gr(1),GX(1)],[gr(2),GX(2)],'-k')
    for ix = ixsg:ixeg
        for iy = iysg:iyeg
            plot(ix,iy,'ko')
        end
    end
    plot([0 0],[0,ny-1],'k')
    plot([nx-1,nx-1],[0,ny-1],'k')
    title(i)
    axis equal
    grid on
    drawnow
    saveas(gcf,strcat(num2str(i),'.jpg'));
    
    clf
end

