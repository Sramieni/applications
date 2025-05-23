function [ xg, Kg ] = newtonBPproto( R, n, coords , xg, Kg)
%newtonBPproto Prototype Newton bifurcation finder

assert(isvector(xg));
assert(~isnan(xg));

import orientationSpace.diffusion.*;

% if(nargout == 0)

Kplot = 8:-0.01:1;
if(isempty(get(groot,'CurrentFigure')) || ~ishold)
out = interpft_extrema(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kplot)); out = orientationSpace.diffusion.alignExtrema(out);
[outdmax,outdmin,~,~,outdother] = interpft_extrema(R.getDerivativeResponseAtPoint(coords.r(n),coords.c(n),1,Kplot));
outdmax = orientationSpace.diffusion.alignExtrema(outdmax);
outdmin = orientationSpace.diffusion.alignExtrema(outdmin);
figure; plot(Kplot,out.');
hold on;
plot(Kplot,outdmax,'k');
plot(Kplot,outdmin,'k');
plot(Kplot,outdother,'ko');
title(sprintf('Local maxima trace for r=%d, c=%d, m=%d',coords.r(n),coords.c(n),coords.m(n)));
xlabel('K');
ylabel('2\theta (Orientation, radians)');
out(:,1:3);
grid on
end
plot(Kg,xg,'ko');

% end

dnK_dmn(:,:,1) = 1;

% xg = maxima2(n)*2; Kg = 8;

last = Inf;

D = 2*pi.^2;

maxKdelta = 4;
jump_threshold = 0.5;

while(abs(dnK_dmn(:,:,1)) > 1e-10)

xold = xg; Kold = Kg;
% [~,dtn_dnm,dnK_dmn] = orientationMaximaTimeDerivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),Kg,3,xg,2*pi,false);
%  || abs(dnK_dmn(:,:,1)./dnK_dmn(:,:,2)) < 0.2 && abs(dnK_dmn(:,:,1)) < 10
% if(abs(dnK_dmn(:,:,1)./dnK_dmn(:,:,2)) < 0.1)
 xgd = interpft1_derivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),xg,2:6);
% d_p2rho_pm2_K = -D.*xgd(:,:,3-1).^2./xgd(:,:,2-1) + D.*xgd(:,:,4-1);
% d_p2rho_pm2_K = d_p2rho_pm2_K.*-4./(2*Kg+1).^3;

dnm_dtn(:,:,1) = -D.*xgd(:,:,3-1)./xgd(:,:,2-1);
dnm_dtn(:,:,2) =           xgd(:,:,3-1) .* dnm_dtn(:,:,1).^2 ...
            +  2 .* D   .* xgd(:,:,4-1) .* dnm_dtn(:,:,1)    ...
            +       D^2 .* xgd(:,:,5-1);
dnm_dtn(:,:,2) = -dnm_dtn(:,:,2)./xgd(:,:,2-1);

 d_p2rho_pm2_t  =                  xgd(:,:,3-1) .* dnm_dtn(:,:,1)    ...
                    +       D   .* xgd(:,:,4-1);
d2_p2rho_pm2_t2 =                  xgd(:,:,3-1) .* dnm_dtn(:,:,2)    ...
                    +              xgd(:,:,4-1) .* dnm_dtn(:,:,1).^2 ...
                    +  2 .* D   .* xgd(:,:,5-1) .* dnm_dtn(:,:,1)    ...
                    +       D^2 .* xgd(:,:,6-1);
                
dnt_dKn(:,:,1) = -4./(2*Kg+1).^3;
dnt_dKn(:,:,2) = 24./(2*Kg+1).^4;

 d_p2rho_pm2_K  = d_p2rho_pm2_t   .* dnt_dKn(:,:,1);
d2_p2rho_pm2_K2 = d2_p2rho_pm2_t2 .* dnt_dKn(:,:,1).^2 ... 
               +  d_p2rho_pm2_t   .* dnt_dKn(:,:,2);

K_jump_estimate = xgd(:,:,2-1) ./ d_p2rho_pm2_K ;
K_jump_estimate2 = 2*xgd(:,:,2-1).*d_p2rho_pm2_K./(2*d_p2rho_pm2_K.^2-xgd(:,:,2-1).*d2_p2rho_pm2_K2);
K_jump_estimate3 = d_p2rho_pm2_K./d2_p2rho_pm2_K2;

% if(K_jump_estimate > 0 && K_jump_estimate2 > 0)
    K_jump_estimate = min(K_jump_estimate,K_jump_estimate2);
% end
if(d2_p2rho_pm2_K2 < 0 && K_jump_estimate > jump_threshold)
    K_jump_estimate = jump_threshold;
end
if( K_jump_estimate < jump_threshold && K_jump_estimate >= 0 && ~(K_jump_estimate3 >= 0))
% if(true)
    disp('Theta Derivative based jump');
    [~,dtn_dnm,dnK_dmn] = orientationMaximaTimeDerivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),Kg,4,xg,2*pi,false);
    xg_newton = dnK_dmn(:,:,1)./dnK_dmn(:,:,2);
    xg_halley = 2*dnK_dmn(:,:,1).*dnK_dmn(:,:,2)./(2*dnK_dmn(:,:,2).^2-dnK_dmn(:,:,1).*dnK_dmn(:,:,3));
    if(abs(xg_newton) < abs(xg_halley))
        xg = xold - xg_newton;
        plot([Kold Kold],[xold xg],'b--');
    else
        xg = xold - xg_halley;
        plot([Kold Kold],[xold xg],'r--');
    end
%      2*v.*vd./(2*vd.^2-v.*vdd);
    Kgpd = (xg-xold).*dnK_dmn(:,:,1)+(xg-xold).^2.*dnK_dmn(:,:,2)/2+(xg-xold).^3.*dnK_dmn(:,:,3)/6+(xg-xold).^4.*dnK_dmn(:,:,4)/24;
    Kg = NaN;
    if(Kgpd < 0 && abs(Kgpd) < 1)
        test = linspace(xold,xg,100);
        Ktest = Kold+(test-xold)*dnK_dmn(:,:,1)+(test-xold).^2*dnK_dmn(:,:,2)/2+(test-xold).^3*dnK_dmn(:,:,3)/6+(test-xold).^4.*dnK_dmn(:,:,4)/24;
        plot(Ktest,test,'g:');
        Kg = halleyK(xg,Kold+Kgpd,R,coords.r(n),coords.c(n));
        plot([Kold Kold+Kgpd],[xg xg],'y:');
        plot([Kold+Kgpd Kg],[xg xg],'m:');
    end
    if(isnan(Kg))
        Kg = halleyK(xg,Kold,R,coords.r(n),coords.c(n));
        plot([Kold Kg],[xg xg],'c:');
    end
    if(abs(Kg - Kold) > jump_threshold)
        Kg = NaN;
    end
    if(Kg - Kold > 1e-10)
        Kg = NaN;
    end
    if(isnan(Kg))
        jump_threshold = jump_threshold/2;
        fprintf('jump threshold now %d\n',jump_threshold);
        Kg = Kold;
        xg = xold;
    else
        jump_threshold = min(jump_threshold*2,0.5);
    end
    plot(Kg,xg,'.');
    if(abs(xgd(:,:,2-1)) > abs(last))
        break;
    else
        last = xgd(:,:,2-1);
    end
%     if(abs(last) < abs(dnK_dmn(:,:,1)))
%         dnK_dmn(:,:,1)
% 
%         break;
%     else
%         last = dnK_dmn(:,:,1);
%     end
else
%     dnm_dKn(:,:,1) = 1./dnK_dmn(:,:,1);
%     dnm_dKn(:,:,2) = -dnK_dmn(:,:,2).*dnm_dKn(:,:,1).^3;
%     dnm_dKn(:,:,3) = (-3*dnK_dmn(:,:,2).*dnm_dKn(:,:,2)-dnK_dmn(:,:,3).*dnK_dmn(:,:,1).^(-2)).*dnK_dmn(:,:,1).^(-2);
%     Kg = Kg - 0.5;
%     xg = halleyft( R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg), xg - 0.5./dnK_dmn(:,:,1),false,1);
    disp('K jump');
    if(K_jump_estimate < 0)
%         K_jump_estimate = maxKdelta;
        K_jump_estimate = d_p2rho_pm2_K./d2_p2rho_pm2_K2;
        if(K_jump_estimate < 0)
            K_jump_estimate = 0.1;
        else
            K_jump_estimate = min(K_jump_estimate*2,jump_threshold);
        end
    end
    xg_change = -Inf;
    xg = NaN;
    while(isnan(xg) || min(xg_change,2*pi-xg_change) > pi/12)
        Kg = Kold - min(K_jump_estimate,maxKdelta);
        Kg = max(Kg,0);
        xg = halleyft( R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg), xold,false,1);
        xg(interpft1_derivatives(R.getResponseAtOrderFTatPoint(coords.r(n),coords.c(n),Kg),xg,2,2*pi) > 0) = NaN;
        xg_change = abs(xg - xold);
        if(isnan(xg) || min(xg_change,2*pi-xg_change) > pi/12)
            maxKdelta = maxKdelta/2;
            fprintf('maxKdelta now %d\n',maxKdelta)
            Kg = Kold;
%             xg = xold;
        else
            plot([Kold Kg],[xold xold],'b-.');
            plot([Kg Kg],[xold xg],'r-.');
        end
    end
%     if(Kg < 1)
%         Kg = Kold;
%         xg = xold;
%         break;
%     end
    plot(Kg,xg,'.');
    last = -Inf;
    jump_threshold = min(jump_threshold*2,0.5);
end

if(Kg <= 0)
    dnK_dmn(:,:,1);

    break;
end



fprintf('xg = %d , Kg = %d, xgd2 = %d\n',xg,Kg,xgd(:,:,2-1));
dnK_dmn(:,:,1);

end

fprintf('xg = %d , Kg = %d, xgd2 = %d\n',xg,Kg,xgd(:,:,2-1));
plot(Kg,xg,'s');



end

