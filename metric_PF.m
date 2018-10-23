function A = metric_PF(x,m,kernel,T)
% Core implementation of our metric:
% Isao Ishikawa, Keisuke Fujii, Masahiro Ikeda, Yuka Hashimoto,and Yoshinobu Kawahara, 
% "Metric on nonlinear dynamical systems with Perron-Frobenius operators," 
% in Advances in Neural Information Processing Systems (Proc. of NIPS), vol. 31, 2018 (to appear).
% arXiv preprint: <https://arxiv.org/abs/1710.04340>

% x: input data
% m: order of outer product
% kernel: type of kernel
% T: time length

if nargin == 3
    [T, N] = size(x{1}) ; % time, dimension
elseif nargin == 4
    N  = size(x{1},2) ;
end

for k = 1:length(x) % 
    for l =  1:length(x) %
        if k <= l 
            sig = MedianDist([x{k};x{l}]) ;
            if m == 1
                % compute kernel
                for t = 1:T
                    for s = 1:N
                        if strcmp(kernel,'Gauss')
                            ker(k,k,t,s) = exp(-abs(x{k}(t,s)-x{k}(t,s)).^2/sig) ;
                            ker(l,l,t,s) = exp(-abs(x{l}(t,s)-x{l}(t,s)).^2/sig) ;
                            ker(k,l,t,s) = exp(-abs(x{k}(t,s)-x{l}(t,s)).^2/sig) ;
                        elseif strcmp(kernel,'Hardy')
                            ker(k,k,t,s) = 1/(1-x{k}(t,s)*conj(x{k}(t,s))) ;
                            ker(l,l,t,s) = 1/(1-x{l}(t,s)*conj(x{l}(t,s))) ;
                            ker(k,l,t,s) = 1/(1-x{k}(t,s)*conj(x{l}(t,s))) ;
                        end
                    end
                end
                % compute dT
                dT(k,k) = sum(sum(ker(k,k,:,:),4),3) ;
                dT(l,l) = sum(sum(ker(l,l,:,:),4),3) ;
                dT(k,l) = sum(sum(ker(k,l,:,:),4),3) ;
                
            elseif m >= 2
                % compute kernel
                for t1 = 1:T
                    for t2 = 1:T
                        for s1 = 1:N
                            for s2 = 1:N
                                if strcmp(kernel,'Gauss')
                                    ker(k,k,t1,t2,s1,s2) = exp(-abs(x{k}(t1,s1)-x{k}(t2,s2)).^2/sig) ;
                                    ker(l,l,t1,t2,s1,s2) = exp(-abs(x{l}(t1,s1)-x{l}(t2,s2)).^2/sig) ;
                                    ker(k,l,t1,t2,s1,s2) = exp(-abs(x{k}(t1,s1)-x{l}(t2,s2)).^2/sig) ;
                                elseif strcmp(kernel,'Hardy')
                                    ker(k,k,t1,t2,s1,s2) = 1/(1-x{k}(t1,s1)*conj(x{k}(t2,s2))) ;
                                    ker(l,l,t1,t2,s1,s2) = 1/(1-x{l}(t1,s1)*conj(x{l}(t2,s2))) ;
                                    ker(k,l,t1,t2,s1,s2) = 1/(1-x{k}(t1,s1)*conj(x{l}(t2,s2))) ;
                                end
                            end
                        end
                    end
                end
                % combination
                if N == 1 % dim = 1 
                    C = [1 1] ;
                else ; C = nchoosek(1:N,2) ; 
                end
                
                % compute IT
                if m == 2
                    for c = 1:size(C,1)
                        for t1 = 1:T
                            for t2 = 1:T
                                matkk = [ker(k,k,t1,t1,C(c,1),C(c,1)) ker(k,k,t1,t2,C(c,1),C(c,2)) ;...
                                    ker(k,k,t2,t1,C(c,2),C(c,1)) ker(k,k,t2,t2,C(c,2),C(c,2))] ;
                                matll = [ker(l,l,t1,t1,C(c,1),C(c,1)) ker(l,l,t1,t2,C(c,1),C(c,2)) ;...
                                    ker(l,l,t2,t1,C(c,2),C(c,1)) ker(l,l,t2,t2,C(c,2),C(c,2))];
                                matkl = [ker(k,l,t1,t1,C(c,1),C(c,1)) ker(k,l,t1,t2,C(c,1),C(c,2)) ;...
                                    ker(k,l,t2,t1,C(c,2),C(c,1)) ker(k,l,t2,t2,C(c,2),C(c,2))];
                                det2(k,k,t1,t2,c) = det(matkk) ;
                                det2(l,l,t1,t2,c) = det(matll) ;
                                det2(k,l,t1,t2,c) = det(matkl) ;
                            end
                        end
                    end
                    dT(k,k) = sum(sum(sum(det2(k,k,:,:,:),5),4),3) ;
                    dT(l,l) = sum(sum(sum(det2(l,l,:,:,:),5),4),3) ;
                    dT(k,l) = sum(sum(sum(det2(k,l,:,:,:),5),4),3) ;
                else error('m > 2 is unsupported');
                end
            end
            % compute A
            A(k,l) = abs(abs(dT(k,l))^2/dT(k,k)/dT(l,l));
        else
            A(k,l) = A(l,k);
        end
    end
end

    function sig = MedianDist(X)
        %  Computing the median of the distances from a data matrix.
        %   Date: Feb 17, 2010
        %   Version: 0.9
        %   Author: Kenji Fukumizu
        %   Affiliation:  The Institute of Statistical Mathematics, ROIS
        %   (c) Kenji Fukumizu
        %-------------------------------------------------------
        %
        %   sig = MedianDist(X)
        %
        %   sig: median of the pairwise distances \|X_i - X_j\|
        %   X:   data matrix
        NN=length(X(:,1));
        ab=X*X';
        aa=diag(ab);
        Dx=repmat(aa,1,NN) + repmat(aa',NN,1) - 2*ab;
        Dx=Dx-diag(diag(Dx));
        dx=nonzeros(reshape(Dx,NN*NN,1));
        sig=sqrt(median(dx));
    end
end