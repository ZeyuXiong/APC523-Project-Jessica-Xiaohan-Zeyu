lis=[9.93186398345891e-07
9.898618625880928e-07
0.00020123340599371986
0.018059756907933378
0.0475050508759247
0.06216727556433868
];

sig=[0.01
0.05
0.10
0.25
0.50
1];

semilogy(sig,lis,'--rs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('sigma')
ylabel('Numerical calculation error')
set(gca,'FontSize',20)
% loglog(sig,lis,'--rs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])