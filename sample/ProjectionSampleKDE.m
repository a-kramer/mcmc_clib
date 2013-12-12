m=26;
fid=fopen('ODEmodel11S26P4U_2013-09-12T11h44m.double','r');
Sample=fread(fid,[m+1,inf],'double');
fclose(fid);

N=columns(Sample);
n=1e3;
Theta=Sample([23,26],1:n:N);
logP=Sample(1:m+1,1:n:N);
clear Sample



K=mkde(Theta);
Mean=[mean(Theta(1,:)),mean(Theta(2,:))];
SD=[std(Theta(1,:)),std(Theta(2,:))];
a=2.3;
ax=[Mean(1)-a*SD(1),Mean(1)+a*SD(1),Mean(2)-a*SD(2),Mean(2)+a*SD(2)];
nx=128;
ny=nx;
x=linspace(ax(1),ax(2),nx);
y=linspace(ax(3),ax(4),ny);

X=cell(1,2);
[X{:}]=meshgrid(x,y);

P=zeros(nx,ny);
for i=1:nx
 for j=1:ny
  theta=[X{1}(i,j);X{2}(i,j)];
  P(i,j)=K(theta);
 endfor
endfor

figure(1); 
plot(Theta(1,:),Theta(2,:),'+;;','color',[0.5,0.5,0.5],'markersize',0.7); hold on;
contour(X{:},P,5,'linewidth',3);
axis(ax*1.1,'equal');
xlabel('$\theta_{23}$');
ylabel('$\theta_{26}$');
hold off

set(gcf,'papersize',[4,4]/2.54);
set(gcf,'paperposition',[0,0,4,4]/2.54);
Folder='~/Dropbox/OxfordBioInf_MCMC_CLIB_paper'
Name='ProjectedSample';
print('-depslatex','-mono',sprintf("%s.tex",Name));
system(sprintf("epstopdf %s.eps",Name));
system(sprintf("mv %s.{tex,pdf} %s",Name,Folder));

