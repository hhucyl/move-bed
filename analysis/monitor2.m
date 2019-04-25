clear
clc
load('2.mat');
contour(ffi);
fffi = mean(ffi,1);
figure
temp = fffi(11:end-10);
plot(temp)
kk = find(temp<=0.1);
kk(1)
Vv = mean(V,2);
figure
plot(Vv,1:numel(Vv))
Q = sum(Vv);
hb = kk(1);
hf = numel(Vv)-hb;

Reb = Q/0.01
2.4*(hf/20)^2

figure
plot(pp(:,hb))

