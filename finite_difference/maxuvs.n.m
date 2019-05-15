lis =[8.353653890391739
2.5294012285957264
0.8322326023302252
0.6433415197596879
0.6294549806830885
0.6270603991432508
0.6265640840275644
0.6264281590815105
];
N=[10
20
40
80
160
250
300
320
];


loglog(N,lis,'--rs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
ylim([0.5 10])
xlim([10 1280])
xlabel('Interval number N')
ylabel('Maximum u')
set(gca,'FontSize',20)