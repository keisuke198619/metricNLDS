% MedianDist()
% 
%  Computing the median of the distances from a data matrix. 
%
%
%   Date: Feb 17, 2010
%   Version: 0.9
%   Author: Kenji Fukumizu
%   Affiliation:  The Institute of Statistical Mathematics, ROIS
%   (c) Kenji Fukumizu
%-------------------------------------------------------
%
%   s = MedianDist(X)
%
%   s: median of the pairwise distances \|X_i - X_j\|
%   X: data matrix
%


function s = MedianDist(X)

N=length(X(:,1)); 
ab=X*X'; 
aa=diag(ab);
Dx=repmat(aa,1,N) + repmat(aa',N,1) - 2*ab;  
Dx=Dx-diag(diag(Dx));
dx=nonzeros(reshape(Dx,N*N,1));
s=sqrt(median(dx));
    