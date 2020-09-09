% Copyright by Yoshinobu Kawahara, Sept. 2016
% Please cite: Yoshinobu Kawahara, "Dynamic Mode Decomposition
% with Reproducing Kernels for Koopman Spectral Analysis," in 
% Advances in Neural Information Processing Systems 29, 2016.

function [lam,Ti,M,F,evec,a,H] = kdmd(Gxx,Gxy,eps)

n = size(Gxx,1)+1;

vOne = ones(n-1,n-1)/(n-1); 
cG   = Gxx-vOne*Gxx-Gxx*vOne+vOne*Gxx*vOne; % centered gram matrix

[V,D]  = eig(cG/(n-1)); 
eval   = diag(D); 
[~,IX] = sort(eval,'descend');
IX     = IX(abs(eval)>eps); 
evec   = V(:,IX); 
n_evec = sqrt(sum(evec.^2));
evec   = evec./repmat(n_evec,size(evec,1),1); 

M = evec*pinv(sqrtm((n-1)*D(IX,IX))); % BS^(1/2) 
H = eye(n-1,n-1)-ones(n-1,n-1)/(n-1); 
F = M'*H*Gxy*H*M; %
[Ti,lam,a] = eig(F); % 
lam = diag(lam); %