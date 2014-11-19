clear
clc
clf
bdgE = load('edge_spectrum.OUT');
AKX = load('AKX.OUT');
figure(1)
T = 26.24;
set(gca,'fontsize',16)
for i = 1:length(bdgE(1,:))
    if max(abs(bdgE(:,i)))<2*pi/T
    plot(AKX,bdgE(:,i)/(pi/T))
    hold on
    end
end
hold off
xlabel('k_x/k_F')
ylabel('\epsilon/(\pi/T)')
axis([-0.5 0.5 -1 1])
load('bdgR.OUT')
load('bdgI.OUT')
bdgH  = bdgR+1i*bdgI;
figure(2)
spy(bdgH)