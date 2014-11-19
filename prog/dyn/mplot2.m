clear
clc
clf
bdgE = load('bulk_spectrum.OUT');
curv = load('curvature.OUT');
AKX = load('AKX.OUT');
AKY = load('AKY.OUT');
figure(1)
mesh(AKX,AKY,curv)
figure(2)
T = 26.8685;
set(gca,'fontsize',16)
for i = 1:length(bdgE(1,:))
    temp = reshape(bdgE(:,i),length(AKX),length(AKX));
    
    mesh(AKX, AKY,temp/(pi/T))
    hold on

end
hold off
xlabel('k_x/k_F')
ylabel('k_y/k_F')
zlabel('\epsilon/(\pi/T)')
axis([-0.2 0.2 -0.2 0.2 -2 2])
view([0 0])
