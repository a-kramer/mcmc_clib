function [kde]=mkde(X)
% Copyright (C) 2010 Andrei Kramer <andrei.kramer@ist.uni-stuttgart.de>
%
% Usage: [kde]=mkde(X)
% X: Sample
%    either rows or columns represent sample members
% returns a function: kde(x), where x is a column vector
% 
 [m,n]=size(X);
 if (m>n)
         X=X';
     [m,n]=size(X);
 end%if
 C=cov(X')*n^(-2/(m+4));
 B=sqrtm(C);
 D2=sqrt(2*pi)^m*sqrt(det(C));
 kde=@(x) mean(exp(-0.5*sumsq(B\bsxfun(@minus,x,X),1))/D2);
end%function
