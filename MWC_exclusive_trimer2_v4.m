clear all

global  cdu_i cdf_i el_i kr_i nloops

nloops = 2;
unit = 1e6;    

if nloops==1
    load trimer_1_loop.txt
    data = trimer_1_loop;
    fprintf('Trimer 1 loop:\n'),
    figure(1), title 'Trimer 1 loop'
    cdf_i = 1; cdu_i = 0;
    el_i = 10;
    kr_i = 1e6/.22;
    
elseif nloops==2
    load trimer_2_loops.txt
    load D3T15T7T8_NaP_37C_all_norm.txt
    %data = trimer_2_loops(1:end,:);
    data = D3T15T7T8_NaP_37C_all_norm(1:end,:);
    fprintf('Trimer 2 loops:\n')
    figure(2), title 'Trimer 2 loops'
    cdf_i = 0; cdu_i = 1;
    el_i = 1/50;
    kr_i = 1/(unit*.184e-6);
    
end
    
lig = data(:,1)*unit;
nn = length(lig);

for i = 2:length(data(1,:))

    F = data(:,i);

    if lig(1) == 0, lig(1) = 1e-9*unit; end
    semilogx(lig,F,'ro','MarkerFaceColor','r'), hold on
    %semilogx(lig,F,'o'), hold on
    lig(1) = data(1,1)*unit;

    [fitresult,gof,out] = createFit(lig,F);
    coeffs = coeffvalues(fitresult);
    cdf = coeffs(1); cdu = coeffs(2); el = coeffs(3); kr = coeffs(4);
    cdf_all(i-1) = cdf; cdu_all(i-1) = cdu; el_all(i-1) = el; kr_all(i-1) = kr;
    kr = kr*unit; % the allosteric constant always comes in molar units, and so should the affinity and disociation constants
    Kd1_all(i-1) = (1+el+el^2)/(el*kr*(1+el)); Kd2_all(i-1) = (1+el*(1+el)*(1+kr))/(el^2*kr^2); Kd3_all(i-1) = (1+el*(1+kr+el*(1+kr+kr^2)))/(el^2*kr^3);
    kr = kr/unit;
    
    CI = confint(fitresult);
    CI(:,4) = 1./(CI(:,4));
    CI = abs(CI(2,:) - CI(1,:))/2;
    e_cdf = CI(1); e_cdu = CI(2); e_el = CI(3); e_kr = CI(4);

    R = out.residuals;
    J = out.Jacobian;
    N = out.numobs; p = out.numparam;
    MSE = (R'*R)/(N-p);
    CovB = inv(J'*J)*MSE;
    
    [fitresult_Hill,gof_Hill,out_Hill] = createFit_Hill(lig,F);
    coeffs_Hill = coeffvalues(fitresult_Hill);
    cdf_Hill = coeffs_Hill(1); cdu_Hill = coeffs_Hill(2); kd_Hill = coeffs_Hill(3); n_Hill = coeffs_Hill(4);

    ligsim = 10.^(-9:0.05:-4)*unit;
    Fsim = cdu + cdf*kr*ligsim.*(el+el^2+2*el^2*kr*ligsim+el^2*kr^2*ligsim.^2)./...
        (1+3*el+el^2 + kr*ligsim.*(3*el+3*el^2+3*el^2*kr*ligsim+el^2*kr^2*ligsim.^2));
    FHill = (cdf_Hill*ligsim.^n_Hill+cdu_Hill*kd_Hill^n_Hill)./(kd_Hill^n_Hill+ligsim.^n_Hill);
    semilogx(ligsim,Fsim,'r-',ligsim,FHill,'r--'), hold off
    %semilogx(ligsim,Fsim,'-')

    fprintf('kr = %g ± %g uM\n',1e6/kr/unit,e_kr*1e6/unit)
    fprintf('el = %g ± %g\n',el,e_el)
    fprintf('cdf = %g ± %g    cdu = %g ± %g\n',cdf,e_cdf,cdu,e_cdu)
    fprintf('rmse = %g\n',gof.rmse)
    fprintf('Kd1 = %g\tKd2 = %g\tKd3 = %g uM\n\n',Kd1_all(i-1)*1e6,Kd2_all(i-1)*1e6,Kd3_all(i-1)*1e6)
        

end

cdf = mean(cdf_all); cdu = mean(cdu_all); el = mean(el_all); kr = mean(1e6./kr_all)/unit;
e_cdf = std(cdf_all); e_cdu = std(cdu_all); e_el = std(el_all); e_kr = std(1e6./kr_all)/unit;
Kd1 = mean(Kd1_all); Kd2 = mean(Kd2_all); Kd3 = mean(Kd3_all);
e_Kd1 = std(Kd1_all); e_Kd2 = std(Kd2_all); e_Kd3 = std(Kd3_all);

fprintf('---------------------------------------\n')
fprintf('kr = %.3g ± %.3g uM\n',kr,e_kr)
fprintf('el = %.3g ± %.3g\n',el,e_el)
fprintf('cdf = %.3g ± %.3g    cdu = %.3g ± %.3g\n',cdf,e_cdf,cdu,e_cdu)
fprintf('Kd1 = %.3g ± %.3g\tKd2 = %.3g ± %.3g\tKd3 = %.3g ± %.3g uM\n\n',Kd1*1e6,e_Kd1*1e6,Kd2*1e6,e_Kd2*1e6,Kd3*1e6,e_Kd3*1e6)

title('Trimer 2 loops')
% DynR = 81.^(1/mean(nH));
% patch(mean(Kd)*[1/sqrt(DynR) 1/sqrt(DynR) sqrt(DynR) sqrt(DynR)],[-.05 1.05 1.05 -.05],'r','FaceAlpha',0.1,'EdgeColor','none')
xlabel('Doxorubicin / M'), ylabel('Occupancy'), axis([unit*1e-9 unit*1e-4 -0.05 1.05])
set(gca,'TickDir','out','FontName','Arial','FontSize',22,'YMinorTick','on')
h = gca; h.YAxis.MinorTickValues = 0.5;
xticks([1e-9 1e-8 1e-7 1e-6 1e-5 1e-4]*unit)
yticks([0 1]), ytickformat('%.1f')
pos = get(gcf,'Position'); set(gcf,'position',[pos(1) pos(2) 450 450])
box off
hold off

figure(4), yyaxis('left')
loglog(ligsim,Fsim./(1-Fsim),'.')
ylabel 'Y / (1 - Y)'
yyaxis('right')
D = gradient(log(Fsim./(1-Fsim)))./gradient(log(ligsim));
plot(ligsim,abs(D),'.')
ylim([0 3]), ylabel 'Hill coefficient'
xlabel 'Ligand'


function [fitresult, gof, out] = createFit(lig, F)

%CREATEFIT(C,CDAVE)
%  Data for 'untitled fit 1' fit:
%      X Input : C
%      Y Output: CDave
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

global  cdf_i cdu_i el_i kr_i  nloops

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( lig, F );

% Set up fittype and options.
MWC_nonexclusive = fittype( 'cdu + cdf*kr*lig.*(el+el^2+2*el^2*kr*lig+el^2*kr^2*lig.^2)./(1+3*el+el^2 + kr*lig.*(3*el+3*el^2+3*el^2*kr*lig+el^2*kr^2*lig.^2))',...
   'independent', 'lig', 'dependent', 'F' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [cdf_i cdu_i el_i kr_i];
opts.MaxFunEvals = 1e6; opts.MaxIter = 1e6;
opts.TolFun = 1e-9; opts.TolX = 1e-9;
opts.Algorithm = 'Trust-Region';
if nloops == 1
    %opts.Lower = [-Inf -Inf 0 1e6/.22 ]; %trimer 1
    %opts.Upper = [Inf Inf Inf 1e6/.22 ]; %trimer 1
    %opts.Lower = [-Inf -Inf 0 1e6/.22 ]; %trimer 1
    %opts.Upper = [Inf Inf Inf 1e6/.22 ]; %trimer 1
    opts.Lower = [-Inf -Inf 0 0]; %trimer 1
    opts.Upper = [Inf Inf Inf Inf]; %trimer 1
elseif nloops == 2
    opts.Lower = [-Inf -Inf 0 0]; %trimer 2
    opts.Upper = [Inf Inf Inf Inf]; %trimer 2
    %opts.Lower = [-Inf -Inf 1/29.403 0]; %trimer 2
    %opts.Upper = [Inf Inf 1/29.403 Inf]; %trimer 2
end

% Fit model to data.
[fitresult, gof, out] = fit( xData, yData, MWC_nonexclusive, opts );

end

function [fitresult, gof, out] = createFit_Hill(lig, F)

%CREATEFIT(C,CDAVE)
%  Data for 'untitled fit 1' fit:
%      X Input : C
%      Y Output: CDave
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( lig, F );

% Set up fittype and options.
Hill = fittype( '(cdf_Hill*lig.^n_Hill+cdu_Hill*kd_Hill^n_Hill)./(kd_Hill^n_Hill+lig.^n_Hill)',...
   'independent', 'lig', 'dependent', 'F' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 0 400e-7 2];
opts.Lower = [-Inf -Inf 0 0]; %trimer 1
opts.Upper = [Inf Inf Inf Inf]; %trimer 1
opts.MaxFunEvals = 1e6; opts.MaxIter = 1e6;
opts.TolFun = 1e-9; opts.TolX = 1e-9;
opts.Algorithm = 'Trust-Region';

% Fit model to data.
[fitresult, gof, out] = fit( xData, yData, Hill, opts );

end