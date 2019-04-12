clear
clc
prefix = {'/home/pzhang/chen/move-bed/'};
fname = {'2.out'};
nul = 0.05;
nu = 1e-6;
gl = 9.8e-4;
rationu = nu/nul;
ratiol = 4e-3/(4*10);
ratiot = ratiol^2/rationu;
g = gl*ratiol/ratiot^2;
name = strcat(prefix,fname);
data = importdata(char(name));
plot(data.data(:,1).*ratiot,data.data(:,2).*ratiol*100)
hold on
plot(data.data(:,1).*ratiot,data.data(:,3).*ratiol*100)
figure
plot(data.data(:,1).*ratiot,data.data(:,4).*ratiol*100/ratiot)
hold on
plot(data.data(:,1).*ratiot,data.data(:,5).*ratiol*100/ratiot)