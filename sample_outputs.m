T=[0.1,8,24,48,72];

ref_u=[0  0  0  0];

u=[
1  0  0  0;
0  1  0  0;
0  0  1  0;
0  0  0  1
];

C=[
0   0   0   1   1   0   0   0   0   0   0;
0   0   1   1   1   0   0   0   0   0   0;
0   0   0   0   0   1   1   0   0   0   0;
0   0   0   0   0   1   0   0   0   0   0;
0   0   0   0   0   0   0   1   1   1   1;
0   0   0   0   0   0   0   1   0   0   0;
];

ny=rows(C);
nu=rows(u);

data=zeros(length(T),ny,nu);
sd_data=data;

reference_data=[
  1195    2548    3690    4488   10961    5753
  1195    2548    3690    4488   10961    5753
  1195    2548    3690    4488   10961    5753
  1195    2548    3690    4488   10961    5753
  1195    2548    3690    4488   10961    5753
];
sd_reference_data=[
  742   742   742   742   742   742
  742   742   742   742   742   742
  742   742   742   742   742   742
  742   742   742   742   742   742
  742   742   742   742   742   742
];


# Data Block for Input line 1
data(:,:,1)=[
        0        0        0        0        0        0
     1359        0     2844     3436        0        0
     4223     3609     3309     2818    15382     5450
     6509     3534     3371     2977        0        0
     4265     3151     2885     3091        0        0
];

# Data Block for Input line 2
data(:,:,2)=[
        0     2781        0        0        0        0
        0        0        0        0        0        0
        0     2974        0        0    16723     6520
        0     3949        0        0        0        0
        0        0        0        0        0        0
];

# Data Block for Input line 3
data(:,:,3)=[
        0        0     8745        0        0        0
        0        0        0        0        0        0
        0        0    15643        0        0        0
        0        0    58271        0        0        0
        0        0        0        0        0        0
]
# Data Block for Input line 4
data(:,:,4)=[
        0        0        0        0    30743        0
        0        0        0        0        0        0
        0        0        0        0    42680        0
        0        0        0        0   119400        0
        0        0        0        0        0        0
];

sd_data=cat(3,
[    inf     inf     inf     inf     inf     inf
    1174     inf    1174    1174     inf     inf
    1174    1174    1174    1174    1174    1174
    1174    1174    1174    1174     inf     inf
    1174    1174    1174    1174     inf     inf
],
[    inf     934     inf     inf     inf     inf
     inf     inf     inf     inf     inf     inf
     inf     934     inf     inf     934     934
     inf     934     inf     inf     inf     inf
     inf     inf     inf     inf     inf     inf
],
[    inf     inf    3695     inf     inf     inf
     inf     inf     inf     inf     inf     inf
     inf     inf    3695     inf     inf     inf
     inf     inf    3695     inf     inf     inf
     inf     inf     inf     inf     inf     inf
],
[    inf     inf     inf     inf   16555     inf
     inf     inf     inf     inf     inf     inf
     inf     inf     inf     inf   16555     inf
     inf     inf     inf     inf   16555     inf
     inf     inf     inf     inf     inf     inf
]);


data=bsxfun(@rdivide,data,reference_data);
sd_data=bsxfun(@times,data,sd_reference_data./reference_data) + bsxfun(@rdivide,sd_data,reference_data);
clear reference_data;

nx=11;
x0=ones(nx,1)*1000;
x0(5)=0;
nt=256;
t=linspace(-T(end),T(end)*1.1,nt);

N=columns(Theta);

figure(10); clf;

l=ceil(N/2);
[sP,I]=sort(logP);
sP=exp(sP/mean(logP));

NC=128;
CMAP=jet(NC);
Ymax=permute(max(data),[3,2,1])';
Ymin=permute(min(data),[3,2,1])';
Ymax(isna(Ymax))=0;
Ymin(isna(Ymin))=0;


h_err=zeros(ny,nu);

ij=(data==0);
data(ij)=NA;
sd_data(ij)=NA;



for k=1:nu % experiments
 for j=l:N % parameters
  theta=Theta(:,I(j));
  rho=exp(theta);
  f=@(x,t) ODEmodel11S26P4U_vf(t,x,[rho',u(k,:)]);
  Jf=@(x,t) ODEmodel11S26P4U_jac(t,x,[rho',u(k,:)]);
  X=abs(lsode({f,Jf},x0,t));
  ref_f=@(x,t) ODEmodel11S26P4U_vf(t,x,[rho',ref_u]);
  ref_Jf=@(x,t) ODEmodel11S26P4U_jac(t,x,[rho',ref_u]);
  ref_X=abs(lsode({ref_f,ref_Jf},x0,t));

  Y=abs(C*X')./(abs(C*ref_X')+1e-2);

  # select color index according to sorted posterior values
  ci=floor(1+(sP(j)-sP(l))/(sP(N)-sP(l))*(NC-1));

  for i=1:ny
   subplot(nu,ny,(k-1)*ny+i);
   plot(t,Y(i,:),"color",CMAP(ci,:),"linewidth",2); hold on;
   if (max(Y(i,:))>Ymax(i,k))
     Ymax(i,k)=max(Y(i,:));
   endif
   if (min(Y(i,:))<Ymin(i,k))
     Ymin(i,k)=min(Y(i,:));
   endif
   axis([-T(end),T(end),Ymin(i,k),Ymax(i,k)]*1.1);
  endfor
 endfor 
 for i=1:ny
  subplot(nu,ny,(k-1)*ny+i);
  if (any(isfinite(data(:,i,k))))
       h_err(i,k)=errorbar(T,data(:,i,k),sd_data(:,i,k));
   set(h_err(i,k),"linewidth",6,"marker",".","linestyle","none","color",[1,1,1]);
       h_err(i,k)=errorbar(T,data(:,i,k),sd_data(:,i,k));
   set(h_err(i,k),"linewidth",2,"marker",".","linestyle","none","color",[0,0,0]);
  endif
   set(gca,"xtick",[0.1,8,24,48,72]);
   set(gca,"xticklabel",{".1","8","24","48","72"});
 endfor
endfor
hold off;

set(gcf,"papersize",[8,6]);
set(gcf,"paperposition",[0,0,8,6]);
set(gcf,"paperorientation","landscape");
print -depslatex SampleOutputsTauUpper.tex
system("epstopdf SampleOutputsTauUpper.eps && rm SampleOutputsTauUpper.eps;")
