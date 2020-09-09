function kpa = kernel_KDMD(dat,MedDistLen)
% if rapidly converging to zero, MedDistLen (e.g, 10 timestep) should be set (otherwise,[])

ns = length(dat) ;
OriginMDL = MedDistLen ;
for n = 1:ns
    data = dat{n} ;
    
    X = data(1:end-1,:) ;
    Y = data(2:end,:) ;
    [tau,dim] = size(X);
    
    % MedDistLen is time length if not converging to zero
    if isempty(OriginMDL)
        MedDistLen = size(X,1);
    elseif OriginMDL > size(X,1)
        MedDistLen = size(X,1);
    end
    
    sig = 2*MedianDist(X(1:MedDistLen,:)); %
    
    % construct gram matrices
    Gxx = gram_gaussian2(X,X,sig);
    Gxy  = gram_gaussian2(X,Y,sig);
    
    % perform kdmd
    [lamn{n},Tin{n},Mn{n},Fn{n},Vn{n},an{n},H] = kdmd(Gxx,Gxy,1e-5);
    
    pn(n) = length(lamn{n});
    Xn{n} = X ;
    Hn{n} = H ;
    Gxxn{n} = Gxx ;
    Gxyn{n} = Gxy ;
    lx(n) = length(Gxx) ;
    sign{n} = sig ;
    disp(['kdmd ',' n=',num2str(n)]);
end
kpa{1} = eye(ns) ;
kpa{2} = eye(ns) ;

for n = 1:ns
    if n >= 2
        Ti2 = Tin{n} ;
        M2 = Mn{n} ;
        H2 = Hn{n} ;
        Gxx2 = Gxxn{n} ;
        X2 = Xn{n} ;
        V2 = Vn{n} ;
        [tau2,dim] = size(X2);
        sig2 = sign{n} ;
        lam2 = lamn{n} ;
        p2 = length(lam2);
        a2 = an{n} ;
        BB = Ti2'*M2'*H2*Gxx2*H2*M2*Ti2;
        for j = 1:n-1
            p1 = length(lamn{j});
            Ti1 = Tin{j} ;
            M1 = Mn{j} ;
            H1 = Hn{j} ;
            X1 = Xn{j} ;
            V1 = Vn{j} ;
            [tau1,dim] = size(X1);
            sig1 = sign{j} ;
            Gxx1 = Gxxn{j} ;
            lam1 = lamn{j} ;
            a1 = an{j} ;
            
            TMP   = gram_gaussian2([X1;X2],[X1;X2],(sig1+sig2)/2.0);
            Gxx12 = TMP(1:tau1,tau1+1:tau1+tau2);
            
            AB = Ti1'*M1'*H1*Gxx12*H2*M2*Ti2;
            AA = Ti1'*M1'*H1*Gxx1*H1*M1*Ti1;
            
            % Lagrange formulation
            % see Section 3.1 in Wolf and Shashua, JMLR, 2003,
            % Learning over Sets using Kernel Principal Angles.
            [V,D] = eig([zeros(p2,p2) AB'; AB zeros(p1,p1)],[BB zeros(p2,p1); zeros(p1,p2) AA]);
            lams = sort(abs(diag(D)),'descend');
            lams = lams(1:min(p1,p2));
            if max(lams)~=0
                lams = lams./max(lams) ;
            else
                disp(['lambdas in kpa1 between ',num2str(n),' and ',num2str(j) ' are zero']) ;
            end
            kpa{1}(n,j) = mean(lams.^2) ;
            kpa{1}(j,n) = kpa{1}(n,j);

            % eigendecomposition approach
            % see Section 3.2 in Wolf and Shashua, JMLR, 2003,
            % Learning over Sets using Kernel Principal Angles.
            [VA,DA] = eig(AA);
            [VB,DB] = eig(BB);
            QAQB = diag(1./sqrt(diag(DA)))*VA'*AB*VB*diag(1./sqrt(diag(DB)));
            [U,S,V] = svd(QAQB);
            lams = sort(abs(diag(S)),'descend');
            lams = lams(1:min(p1,p2));
            if max(lams)~=0
                lams = lams./max(lams) ;
            else
                disp(['lambdas in kpa2 between ',num2str(n),' and ',num2str(j) ' are zero']) ;
            end
            kpa{2}(n,j) = mean(lams.^2) ;
            kpa{2}(j,n) = kpa{2}(n,j);
        end
    end
    disp(['KPA',' n=',num2str(n)]);
end
