clear all;
N=319;
sigma=0.1;

dx=1/(N+1);
x=linspace(dx,1-dx,N);
y=x;
[Y,X]=meshgrid(y,x);

ii=1:N;
jj=1:N;
kx=ii*pi;
ky=jj*pi;
[KY,KX]=meshgrid(ky,kx);
K2=-KX.^2-KY.^2;

DST_for=@(u)dst(dst(u)')';
DST_back=@(u)idst(idst(u)')';

delta=-exp(-((X-0.5).^2+(Y-0.5).^2)/2/sigma^2)/2/pi/sigma/sigma;
u=DST_back(DST_for(delta)./K2);

MX=max(max(u));
% disp(MX)
%plot and calculate maimum value

%figure;surf(X,Y,u);
% shading interp;colormap jet;
% saveas(gcf,'N=1280.png')
% MX=max(max(u));
% disp(MX)

%calculate relationship of error and sigma
% sum=0;
% for i=1:N
%     sum = sum + u(0.25*(N+1),i)*dx;
% end
% error=abs(0.068184116-sum);

%calculate convergence error 
%sigma=0.01, convergent value =[0.625450068823669]
%sigama=0.05, [0.369300069460278]
%sigma=0.1 [0.258982277780559]
%sigma=0.25 [0.115846970798426]
%sigma=0.5 [[0.0406704993906991]] [0.0406704997520970]
%converg_error1 = abs(MX-0.625450068823669)/0.625450068823669
%converg_error2 = abs(MX-0.369300069460278)/0.369300069460278
converg_error3 = abs(MX-0.258982277780559)/0.258982277780559
%converg_error4 = abs(MX-0.115846970798426)/0.115846970798426
%converg_error5 = abs(MX-0.0406704997520970)/0.0406704997520970
    
