clear all

global  cdu_i cdf_i el_i kr_i nloops

nloops = 1;
unit = 1e9;    

if nloops==1
    load trimer_1_loop.txt
    data = trimer_1_loop;
    fprintf('Trimer 1 loop:\n'),
    figure(1), title 'Trimer 1 loop'
    cdf_i = 1; cdu_i = 0;
    el_i = 10;
    kr_i = 1e6/(unit*.22);
    
elseif nloops==2
    load trimer_2_loops.txt
    data = trimer_2_loops(1:end,:);
    fprintf('Trimer 2 loops:\n')
    figure(2), title 'Trimer 2 loops'
    cdf_i = 0; cdu_i = 1;
    el_i = 100;
    kr_i = 1/(unit*.22e-6);
    
end
    
lig = data(:,1)*unit;
F = data(:,2);
nn = length(lig);

if lig(1) == 0, lig(1) = unit*1e-9; end
semilogx(lig,F,'ko','MarkerFaceColor','k'), hold on
lig(1) = data(1,1)*unit;

[fitresult,gof,out] = createFit(lig,F);
coeffs = coeffvalues(fitresult);
cdf = coeffs(1); cdu = coeffs(2); el = coeffs(3); kr = coeffs(4);
kr = kr*unit; Kd1 = (el+1)/kr; Kd2 = (el+kr+1)/kr^2; Kd3 = (el+kr+1+kr^2)/kr^3; kr = kr/unit;
%Calculation of Kd,i taking into account degenerate states
%This is not correct because Kd,3 is lower than intrinsic Kd.
%Kd1 = (el+1)/(3*kr); Kd2 = (el+3*kr+1)/(3*kr^2); Kd3 = (el+3*kr+1+3*kr^2)/kr^3;

CI = confint(fitresult);
CI(:,4) = CI(:,4)*unit;
CI(:,5) = (CI(:,3)+1)./CI(:,4);
CI(:,6) = (CI(:,3)+CI(:,4)+1)./CI(:,4).^2;
CI(:,7) = (CI(:,3)+CI(:,4).^2+CI(:,4)+1)./CI(:,4).^3;
CI(:,4) = CI(:,4)/unit;
e = abs(CI(2,:) - CI(1,:))/2; %e(4) = 1./(CI(:,4));
e_cdf = e(1); e_cdu = e(2); e_el = e(3);
e_kr = abs(1/CI(2,4)-1/CI(1,4))/2;
e_Kd1 = e(5); e_Kd2 = e(6); e_Kd3 = e(7);

R = out.residuals;
J = out.Jacobian;
N = out.numobs; p = out.numparam;
MSE = (R'*R)/(N-p);
CovB = inv(J'*J)*MSE;

[fitresult_Hill,gof_Hill,out_Hill] = createFit_Hill(lig,F);
coeffs_Hill = coeffvalues(fitresult_Hill);
cdf_Hill = coeffs_Hill(1); cdu_Hill = coeffs_Hill(2); kd_Hill = coeffs_Hill(3); n_Hill = coeffs_Hill(4);

ligsim = 10.^(-9:0.05:-4)*unit;
Fsim = cdu + cdf*(kr*ligsim+2*(kr*ligsim).^2+(kr*ligsim).^3)./...
    (1+el + 3*kr*ligsim + 3*kr^2*ligsim.^2 + kr^3*ligsim.^3);
FHill = (cdf_Hill*ligsim.^n_Hill+cdu_Hill*kd_Hill^n_Hill)./(kd_Hill^n_Hill+ligsim.^n_Hill);
semilogx(ligsim,Fsim,'k-',ligsim,FHill,'k--'), hold off

fprintf('kr = %.2f ± %.2f uM\n',1e6/kr/unit,e_kr*1e6/unit)
fprintf('el = %g ± %g\n',el,e_el)
fprintf('cdf = %g ± %g    cdu = %g ± %g\n',cdf,e_cdf,cdu,e_cdu)
fprintf('rmse = %g\n',gof.rmse)
fprintf('Kd1 = %.3g ± %.3g\tKd2 = %.3g ± %.3g\tKd3 = %.3g ± %.3g uM\n\n',Kd1,e_Kd1,Kd2,e_Kd2,Kd3,e_Kd3)
%fprintf('Kd1 = %g\tKd2 = %g\tKd3 = %g uM\n\n',(el+1)/(3*kr),(el+3*kr+1)/(3*kr^2),(el+3*kr+1+3*kr^2)/kr^3)

title('Trimer 1 loop')
% DynR = 81.^(1/mean(nH));
% patch(mean(Kd)*[1/sqrt(DynR) 1/sqrt(DynR) sqrt(DynR) sqrt(DynR)],[-.05 1.05 1.05 -.05],'r','FaceAlpha',0.1,'EdgeColor','none')
xlabel('Doxorubicin / M'), ylabel('Occupancy'), axis([1e-9*unit 1e-4*unit -0.05 1.05])
set(gca,'TickDir','out','FontName','Arial','FontSize',22,'YMinorTick','on')
h = gca; h.YAxis.MinorTickValues = 0.5;
xticks([1e-9 1e-8 1e-7 1e-6 1e-5 1e-4]*unit)
yticks([0 1]), ytickformat('%.1f')
pos = get(gcf,'Position'); set(gcf,'position',[pos(1) pos(2) 450 450])
box off
hold off

figure(3), yyaxis('left')
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
MWC_nonexclusive = fittype( 'cdu + cdf*(kr*lig+2*(kr*lig).^2+(kr*lig).^3)./(1+el + 3*kr*lig + 3*kr^2*lig.^2 + kr^3*lig.^3)',...
   'independent', 'lig', 'dependent', 'F' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [cdf_i cdu_i el_i kr_i];
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
    opts.Lower = [-Inf -Inf 0 1e6/.18]; %trimer 2
    opts.Upper = [Inf Inf Inf 1e6/.18]; %trimer 2
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

% Fit model to data.
[fitresult, gof, out] = fit( xData, yData, Hill, opts );

end