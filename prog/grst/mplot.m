clear
clc
clf
curv = load('curvature.OUT');
bdgE = load('spectrum.OUT');
AKX = load('AKX.OUT');
AKY = load('AKY.OUT');
figure(1)
set(gca,'fontsize',16)
mesh(AKX/pi,AKY/pi,curv)
xlabel('k_x/k_F')
ylabel('k_y/k_F')
zlabel('curvature')
%axis([min(AKX) max(AKX) min(AKY) max(AKY) -1 1 ])
figure(2)
set(gca,'fontsize',16)
for i = 1:length(bdgE(1,:))
    temp = reshape(bdgE(:,i),length(AKX),length(AKX));
    mesh(AKX/pi,AKY/pi,temp/(pi/25))
    hold on
end
hold off
xlabel('k_x/k_F')
ylabel('k_y/k_F')
zlabel('\epsilon')
axis([-1 1 -1 1 -5 5])
view([0 0])