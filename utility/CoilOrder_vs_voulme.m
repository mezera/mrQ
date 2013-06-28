

for nSamples=1:20
    for pOrder=1:3
[~,~, ~,~,VarEx(nSamples,pOrder)]=pdPolyPhantomOrder(nSamples, 32,3,pOrder,500,3);
    end
end

name=(fullfile(mrqRootPath,'PD','figures', 'VarEx_nova32Ch'));
save(name,'VarEx');
figure;plot(VarEx,'-*')
legend('pOrder_1','pOrder_2','pOrder_3')
xlabel('nSamples')
ylabel('VarEx')


