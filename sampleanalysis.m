m=26;
filename="sample/ODEmodel11S26P4U_2013-09-12T11h44m.double"

fid=fopen(filename,"r");
D=fread(fid,[m+1,Inf],"double"); % 26 parameters + log posterior
fclose(fid);
n=columns(D)

figure(1);
plot(D(m+1,:),"-;log posterior;");

#pause;

LP=zeros(3,1);
tau=zeros(2,1);

[LP(1),LP(2),LP(3),tau_int(1),tau_int(2)]=UWerr(D(m+1,:)');
tau_int
skip=ceil(2*sum(tau_int));
selection=skip:skip:n; # remove some points at the beginning of the sample
Theta=D(1:m,selection);  # just in case they are low probability points
logP=D(m+1,selection);   # i.e. the markov chain did not converge to the posterior yet
clear D;

function [FH]=pcplot(Sample,P,varargin)
 [m,n]=size(Sample);
 [s,I]=sort(P);
 if (nargin>2)
  CMAP=varargin{1};
 else
  CMAP=colormap('bone');
 endif
 c=rows(CMAP);
 lsn=linspace(0,1,n);
 lsc=linspace(0,1,c);
 color_order=zeros(n,3);
 for i=1:3
  color_order(:,i)=interp1(lsc,CMAP(:,i),lsn);
 endfor
 %FH=figure(); clf;
 set(gca,"colororder",color_order);
 FH=plot(Sample,"-;;","linewidth",3); 
endfunction


prior_mu=[
  -1.409318
   1.106765
   1.666112
   4.502982
   0.518634
   1.743055
  -2.966790
  -1.382261
   2.607560
   4.343400
  -0.051188
  -0.271598
   1.694231
  -4.693916
  -2.413029
  -2.346162
  -4.907400
  -0.268742
   2.233809
  -1.376299
   6.675061
  -2.260216
  -3.872369
   4.929244
  -0.474706
  -2.392459
];

printf("plotting the sample can take a while (will plot %i lines).\n",length(selection));
fflush(stdout);
theta=mean(Theta');
FH=pcplot(Theta,exp(logP)); hold on
plot(1:m,theta,"-;mean;","color",[1.0,0.4,0.4],"linewidth",3)

lw={4,2};
lc={[0,0,0],[0.7,1.0,0.7]};
for i=1:2
 h_err=errorbar(prior_mu,ones(m,1)/sqrt(0.2));
 set(h_err,"linewidth",lw{i},"color",lc{i},"marker",".","linestyle","none");
endfor
hold off

xlabel('$i$');
ylabel('$\theta_i$');

axis([0.8,m+0.2,min(min(Theta))*0.9,max(max(Theta))*1.1]);
set(gcf,'papersize',[16,10]/2.54);
set(gcf,'paperposition',[0,0,16,10]/2.54);
set(gcf,'paperorientation','landscape');
print -depslatex sample.tex
system("epstopdf sample.eps && rm sample.eps")
