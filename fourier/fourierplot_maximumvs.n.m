lis =[7.16018815261236,2.22772791278283,0.756616042555661,0.625747113445116,0.625450069706174,0.625450068823664,0.625450068823667,0.625450068823667];
N=[10,20,40,80,160,320,640,1280];


loglog(N,lis,'--rs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
ylim([0.5 10])
xlim([10 1280])

hold on
y=0.625450068823667*ones(1, 128);
plot(10:10:1280, y)
xlabel('Interval number N')
ylabel('Maximum u')
set(gca,'FontSize',20)