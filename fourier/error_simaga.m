lis=[1.06970107216497e-08
1.33700936572057e-08
0.000200418154764431
0.0180595260726323
0.0475047862019134
0.0579121068781060
0.0621671710424389];

sig=[0.01
0.05
0.10
0.25
0.50
0.75
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