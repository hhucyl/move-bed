clear
clc
num = [6,7,8,9];
prefix = {'/home/pzhang/chen/move-bed/'};
fname = {'periodic1_'}; 
for i=1:numel(num)
    name = strcat(prefix,fname,num2str(num(i)),'.out');
    data = importdata(char(name));
    x = data.data(:,1);
    y = data.data(:,2);
    plot(data.data(:,1),data.data(:,2),'linewidth',2)
    hold on
    ylim([0 100])
end
plot(max(data.data(:,1)),0.5979*100,'o','markersize',20)