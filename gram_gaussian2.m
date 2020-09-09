function Gxy = gram_gaussian2(X,Y,sig)
% Keisuke Fujii,2017

n = size(X,1);
s2  = 2*sig*sig;
X2 = repmat(sum(X.^2,2),1,n) ; 
Y2 = repmat(sum(Y.^2,2),1,n)'; 
XY = X*Y';
D2 = X2+Y2-2*XY ;
Gxy = exp(-D2./s2);

% for test
% for i = 1:n
%     for j = 1:n
%         Gxy2(i,j) = exp(-(sum((X(i,:)-Y(j,:)).^2,2))/s2);
%     end
% end
% sum(abs(Gxy(:)-Gxy2(:))) % almost the same
