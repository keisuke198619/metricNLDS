%% demo_metric_PF
% a simple demonstration script of our metric using Rotation on the unit disk:
% Isao Ishikawa, Keisuke Fujii, Masahiro Ikeda, Yuka Hashimoto, and Yoshinobu Kawahara, 
% "Metric on nonlinear dynamical systems with Perron-Frobenius operators," 
% in Advances in Neural Information Processing Systems (Proc. of NIPS), vol. 31, 2018.
clear ; close all

% parameters of data -------------------
T = 100 ; % time length
rng(0, 'twister');

x0 = (rand(1,1)-0.5)*2 + (rand(1,1)-0.5)*2*1i ; % initial values
ampx0 = 0.9; % 0.9/0.3/rand to restrict to the unit disk
x0 = repmat(x0/abs(x0)*ampx0,9,1) ; x0 = x0(:) ;

th = [1/3 1/4 pi/3] ; % rotation angle
qdiff = 12 ; % satisfying th1-th2 = p/q (rational number case)
th = repmat(th,3,1) ; th = th(:) ;
amp = [1 0.9 0.3] ; % amplitude
amp = repmat(amp,1,3) ;
for n = 1:length(th)
    alpha{n} = exp(2*pi*1i*th(n))*amp(n) ;
end

describ = {'1: \theta = 1/3, |\alpha| = 1','2: \theta = 1/3, |\alpha| = 0.9','3: \theta = 1/3, |\alpha| = 0.3',...
    '4: \theta = 1/4, |\alpha| = 1','5: \theta = 1/4, |\alpha| = 0.9','6: \theta = 1/4, |\alpha| = 0.3',...
    '7: \theta = \pi/3, |\alpha| = 1','8: \theta = \pi/3, |\alpha| = 0.9','9: \theta = \pi/3, |\alpha| = 0.3',...
    } ;

% create data -------------------

for n = 1:length(alpha)
    for t = 1:T
        x{n}(t,1) = alpha{n}^(t-1)*x0(n) ; %
    end
    x_2{n} = cat(2,x{n}(1:end-1,1),x{n}(2:end,1)) ;
    if 0 % visualize
        figure(1)
        subplot(3,3,n)
        plot(x{n}) ; hold on;
        plot(x0(n),'bo');
        axis equal
        xlim([-1 1]); ylim([-1 1]);
        box off
        title(describ{n}) ;
    end
end

% calculate metric -------------------

if 1 % the Koopman kernel of principal angle (Fujii et al. 2017)
    % https://link.springer.com/chapter/10.1007/978-3-319-71273-4_11
    
    MedDistLen = 10; % if rapidly converging to zero, MedDistLen (e.g, 10 timestep) 
                     % should be set (otherwise,[] in most dataset)

    kpa = kernel_KDMD(x_2,MedDistLen) ;
end

% numerical
if 1
    A1T10 = metric_PF(x,1,'Hardy',10) ;
    A1Tnm = metric_PF(x,1,'Hardy') ; %
    A2Tnm = metric_PF(x_2,2,'Hardy') ;
    A1T_Gauss = metric_PF(x,1,'Gauss') ;
    A2T_Gauss = metric_PF(x_2,2,'Gauss') ;
end

if 1 % analytical (T -> inf)
    Or = 1 ; % order
    for k = 1:length(x)
        for l = 1:length(x)
            if k <= l
                % analytical (T -> inf.) m = 1
                if amp(k) == 1 && amp(l) == 1
                    % r = log(alpha{k}*conj(alpha{l}))/(2*pi*1i) ; % r = th(l)-th(k) = p/q
                    numerator = (1 - abs(x0(k))^2) * (1 - abs(x0(l))^2) ;
                    if (th(k) == th(end) || th(l) == th(end)) && th(k) ~= th(l) % r is an irrational number
                        A1Tan(k,l) = numerator ;
                    else % r is a rational number
                        zconjw = x0(k)*conj(x0(l)) ;
                        if th(k)==th(l)
                            A1Tan(k,l) = 1 ;
                        else A1Tan(k,l) = numerator/(1-zconjw^qdiff) ;
                        end
                    end
                elseif amp(k) == 1 && amp(l) < 1
                    A1Tan(k,l) = 1 - abs(x0(k))^2 ;
                elseif amp(k) < 1 && amp(l) == 1
                    A1Tan(k,l) = 1 - abs(x0(l))^2 ;
                elseif amp(k) < 1 && amp(l) < 1
                    A1Tan(k,l) = 1 ;
                end
            else
                A1Tan(k,l) = A1Tan(l,k);
            end
        end
    end
end 
if 0% plot -------------------
    cmin = 0 ;
    figure(2)
    subplot 171
    h = heatmap(A1T10); 
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_1 (Szego) T=10']) 
    subplot 172
    h = heatmap(A1Tnm); 
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_1 (Szego)   T=',num2str(T)])
    subplot 173
    h = heatmap(A1Tan); 
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_1 (Szego)   T->inf'])
    subplot 174
    h = heatmap(A2Tnm); 
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_2 (Szego)   T=',num2str(T)])
    subplot 175
    h = heatmap(A1T_Gauss);
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_1g (Gauss)   T=',num2str(T)])
    subplot 176
    h = heatmap(A2T_Gauss);
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_2g (Gauss)   T=',num2str(T)])
    subplot 177
    h = heatmap(kpa{2});
    h.ColorLimits = [cmin 1] ;
    h.ColorbarVisible = 'off' ;
    title(['A_kkp (KDMD)   T=',num2str(T)])
    set(gcf,'paperposition',[-6.5316 12.6955 34.0632 4.3089]);
    set(gcf,'OuterPosition',[106.1429 555.8571 1.1269e+03 242.2857]) ;
end