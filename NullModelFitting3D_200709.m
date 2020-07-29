%%Written by Sriram Srikant
%%Dept. of MCB, Harvard University
%%https://github.com/sriramsrikant

%%Curve fitting script to help Hannah Foster with her project. Scientific
%%details below
%%I am going to set this up using the fitting toolbox that MATLAB ships with.

%%Variable index:
%DNALength(experiments) | DNA length used in different experiments
%DNA(samples) | [DNA] in experiment
%vATP(samples) | measured ATPase rate in experiment
%vMaxATP,A_uni & K_Diss are the parameters for the SchemeII fit and
%vMaxLength,K_Length & K_ATP are the parameters for the SchemeIII fit.
%%Specified Variables:
bindingSite = 20; %This is the size of the binding Site you are going to be testing. Important to remember that the system can't use data with DNALength == bindingSite

%%These are the checks before running the fitting.
%Check if the data variables exists.;
if (exist('DNA','var') ~= 1) || (exist('vATP','var') ~= 1) || (exist('DNALength','var') ~= 1)
    disp('Data not present in the appropriate variables');
    return
end
%Check if there is enough data entered into the data set
if (size(DNA,1)<30) || (size(vATP,1)<30 || size(DNALength,1)<30)
    disp('Not enough data entered in DNA, DNALength and vATP. Atleast 3 data-points needed per experiment');
    return
end
%Check if the data entered into DNA, DNALength and vATP are consistent
if (size(DNA,1) ~= size(vATP,1) || size(DNA,1) ~= size(DNALength,1))
    disp('The data-points in DNA, DNALength and vATP are not consistent');
    return
end

%%I am going to also use a model for a non-translocating DNA-binding
%%protein. According to the theory done by Young et al; J. Mol. Biol.
%%(1994),235:1436-1446 they just use a simple Michaelis-Menton kinetics of
%%ATP hydrolysis by the bound complex. This is going to be used for
%%null-model testing against the SchemeIII and SchemeII fit with an F-test
%%since we have different number of parameters in each of the models.
%%SchemeI predictions:
%%v_{Max,ATP} = constant in DNALength
%%K_{ATP} = constant in DNALength

%%vATP = vMaxATP*DNA/(K_ATP + DNA)
%%with parameters vMaxATP, K_ATP. This can be thought as a restriction of
%%the SchemeIIa/b model with A_uni == 0 && K_diss = K_ATP; and a restriction
%%of SchemeIIIa/b with K_Length == 0 && vMaxLength = vMaxATP. I can then
%%use a standard F-test for null hypothesis testing.

%%I am going to finish the next layer of the model to match the theory from
%%Young et al; J. Mol. Biol. (1994),235:1436-1446 in order to test between
%%SchemeII and SchemeIII.

%%SchemeIIa/b predictions:
%%v_{Max,ATP} = [Enzyme]*k = constant in DNALength
%%  where, k = rate of ATP hydrolysis
%%K_{ATP} = A_{uni}/(DNALength-bindingSite) + K_{diss}
%%  where, A_{uni} = 2*k_{i}'/k_{1}' & K_{diss} = k_{-1}/k_{1} & bindingSite is the size of the binding site in bp

%%vATP = vMaxATP*DNA/((A_uni/(DNALength-bindingSite)) + K_diss + DNA)
%%with parameters, vMaxATP, A_uni, K_diss

%%SchemeIIIa/b predictions:
%%v_{Max,ATP} = [Enzyme]*k*DNALength/((2k_{i}^{2}/(k_{i}+k_{-1})*k_{-1}) +
%%DNALength)
%% v_{Max,ATP} = V_{Max,Length}*(DNALength-bindingSite)/(K_{Length}+(DNALength-bindingSite))
%%  where, k = rate of ATP hydrolysis; 
%%K_{ATP} = K_{d}/k_{i} = constant in DNALength
%%  where, K_{d} = dissociation from end-of-lattice; k_{i} = translocation rate & bindingSite is the size of the binding site in bp

%%vATP = vMaxLength*(DNALength-bindingSite)*DNA/(K_Length + (DNALength-bindingSite))*(K_ATP + DNA)
%%with parameters, vMaxLength, K_Length, K_ATP

%%SchemeIVa/b predictions:
%%v_{ATP} = v_{Max,ATP}*DNA/(K_{ATP} + DNALength)
%% v_{Max,ATP} = V_{Max,Length}*(DNALength-bindingSite)/(K_{Length}+(DNALength-bindingSite))
%%  where, k = rate of ATP hydrolysis; 
%%K_{ATP} = K_{ATP,Length}/(DNALength-bindingSite)
%%  where, K_{ATP,Length} = is the constant of the K_{ATP} dependence on DNAlength

%%vATP = vMaxLength*(DNALength-bindingSite)*DNA/((K_Length +
%%(DNALength-bindingSite))*(K_ATPLength/(DNALength-bindingSite) + DNA))
%%with parameters, vMaxLength, K_{Length}, K_{ATP,Length}

%%SchemeV predictions:
%%v_{ATP} = v_{Max,ATP}*DNAeff/(Keff_{ATP} + DNAeff)
%%DNAeff = DNA*DNALength
%%  where, the effective concentration of bindingSites is proportional to DNA*DNALength

%%v_{ATP} = v_{Max,ATP}*DNA*DNALength/(Keff_{ATP} + DNA*DNALength)
%%with parameters, vMaxATP, Keff_{ATP}


%%This is the model fitting with the data that has been input.
%These are the mechanistically detailed kinetic models for the ATPase rate of the translocase dependent on [DNA]

inputPrompt = input('Which scheme do you want to test? 2 or 3 or 4 or 5 or 6:');
%%This is the SchemeII fitting with a pre-set bindingSite, tested against
%%the SchemeI Null model.
if (inputPrompt == 2)
    indices = find(DNALength > bindingSite);    
    expr = ['vMaxATP*DNA/((A_uni/(DNALength-' num2str(bindingSite) ')) + K_diss + DNA)'];
    IIft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','A_uni','K_diss'});
    [IIfitModel,IIfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),IIft,'Start',[1 1 1],'Lower',[0 0 0]);
    confidenceInter = confint(IIfitModel);
    
    expr = 'vMaxATP*DNA/(K_ATP + DNA) + 0*DNALength';
    Ift = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','K_ATP'});
    [IfitModel,IfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),Ift,'Start',[1 1],'Lower',[0 0]);
    confidenceInterNull = confint(IfitModel);
        
    %Plotting of the fit to SchemeII from all the data.
    figure
    hold on
    subplot(2,1,1)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IIfitModel);
    titleString = ['Fitting data to SchemeII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IfitModel);
    titleString = ['Fitting data to Null Model SchemeI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Plotting of the residual histograms.
    figure
    hold on
    subplot(2,1,1)
    plot(IIfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to SchemeII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    plot(IfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to Null Model SchemeI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Calculating the F-statistic of the model compared to the Null model
    %using the standard definition of the F-statistic for a restricted null
    %model. In the statistic below the Null model is always given by 1 and
    %the alternate hypothesis being the model tested as 2.
    %FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2)
    %diff = vATP(:) - IIfitModel(DNALength(:),DNA(:));
    %SS2 = sum(diff(:).^2);
    SS2 = IIfitStatistics.sse;
    df2 = length(unique([DNALength,DNA], 'rows')) - 3;
    %diff = vATP(:) - IfitModel(DNALength(:),DNA(:));
    %SS1 = sum(diff(:).^2);
    SS1 = IfitStatistics.sse;
    df1 = length(unique([DNALength,DNA], 'rows')) - 2;
    FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2);
    pValue = 1-fcdf(FStatistic,(df1-df2),df2);
    
    %I want to display the final fitted parameters of the data and fits I have.
    disp('The models are fit to the SchemeII with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeII:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(IIfitModel.vMaxATP) '(' num2str(confidenceInter(1,1)) ',' num2str(confidenceInter(2,1)) ');'];
    disp(modelString);
    modelString = ['A_{uni} = ' num2str(IIfitModel.A_uni) '(' num2str(confidenceInter(1,2)) ',' num2str(confidenceInter(2,2)) ');'];
    disp(modelString);
    modelString = [' K_{diss} = ' num2str(IIfitModel.K_diss) '(' num2str(confidenceInter(1,3)) ',' num2str(confidenceInter(2,3)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IIfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('The models are fit to the Null Model SchemeI with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeI:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(IfitModel.vMaxATP) '(' num2str(confidenceInterNull(1,1)) ',' num2str(confidenceInterNull(2,1)) ');'];
    disp(modelString);
    modelString = [' K_{ATP} = ' num2str(IfitModel.K_ATP) '(' num2str(confidenceInterNull(1,2)) ',' num2str(confidenceInterNull(2,2)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('By performing an F-test on the models and comparing them, the following is the p-value, i.e. the probability that null model is not rejected.');
    modelString = ['F-Statistic from F-test = ', num2str(FStatistic)];
    disp(modelString);
    modelString = ['p-Value from F-test = ', num2str(pValue)];
    disp(modelString);

%%This is the SchemeIII fitting with a pre-set bindingSite, tested against
%%the SchemeI Null Model.
elseif (inputPrompt == 3)
    indices = find(DNALength > bindingSite);    
    expr = ['vMaxLength*(DNALength-' num2str(bindingSite) ')*DNA/((K_Length + (DNALength-' num2str(bindingSite) '))*(K_ATP + DNA))'];
    IIIft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxLength','K_Length','K_ATP'});
    [IIIfitModel,IIIfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),IIIft,'Start',[1 1 1],'Lower',[0 0 0]);
    confidenceInter = confint(IIIfitModel);
        
    expr = 'vMaxATP*DNA/(K_ATP + DNA) + 0*DNALength';
    Ift = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','K_ATP'});
    [IfitModel,IfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),Ift,'Start',[1 1],'Lower',[0 0]);
    confidenceInterNull = confint(IfitModel);
    
    %Plotting of the fit to SchemeIII from all the data.
    figure
    hold on
    subplot(2,1,1)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IIIfitModel);
    titleString = ['Fitting data to SchemeIII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IfitModel);
    titleString = ['Fitting data to Null Model SchemeI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Plotting of the residual histograms.
    figure
    hold on
    subplot(2,1,1)
    plot(IIIfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to SchemeIII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    plot(IfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to Null Model SchemeI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Calculating the F-statistic of the model compared to the Null model
    %using the standard definition of the F-statistic for a restricted null
    %model. In the statistic below the Null model is always given by 1 and
    %the alternate hypothesis being the model tested.
    %FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2)
    %diff = vATP(:) - IIfitModel(DNALength(:),DNA(:));
    %SS2 = sum(diff(:).^2);
    SS2 = IIIfitStatistics.sse;
    df2 = length(unique([DNALength,DNA], 'rows')) - 3;
    %diff = vATP(:) - IfitModel(DNALength(:),DNA(:));
    %SS1 = sum(diff(:).^2);
    SS1 = IfitStatistics.sse;
    df1 = length(unique([DNALength,DNA], 'rows')) - 2;
    FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2);
    pValue = 1-fcdf(FStatistic,(df1-df2),df2);
    
    %I want to display the final fitted parameters of the data and fits I have.
    disp('The models are fit to the SchemeIII with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeIII:';
    disp(modelString);
    modelString = ['v_{Max,Length} = ' num2str(IIIfitModel.vMaxLength) '(' num2str(confidenceInter(1,1)) ',' num2str(confidenceInter(2,1)) ');'];
    disp(modelString);
    modelString = ['K_{Length} = ' num2str(IIIfitModel.K_Length) '(' num2str(confidenceInter(1,2)) ',' num2str(confidenceInter(2,2)) ');'];
    disp(modelString);
    modelString = ['K_{ATP} = ' num2str(IIIfitModel.K_ATP) '(' num2str(confidenceInter(1,3)) ',' num2str(confidenceInter(2,3)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IIIfitStatistics.rsquare)];
    disp(modelString);
    
    disp('The models are fit to the Null Model SchemeI with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeI:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(IfitModel.vMaxATP) '(' num2str(confidenceInterNull(1,1)) ',' num2str(confidenceInterNull(2,1)) ');'];
    disp(modelString);
    modelString = [' K_{ATP} = ' num2str(IfitModel.K_ATP) '(' num2str(confidenceInterNull(1,2)) ',' num2str(confidenceInterNull(2,2)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('By performing an F-test on the models and comparing them, the following is the p-value, i.e. the probability that null model is not rejected.');
    modelString = ['F-Statistic from F-test = ', num2str(FStatistic)];
    disp(modelString);
    modelString = ['p-Value from F-test = ', num2str(pValue)];
    disp(modelString);

%%This is the SchemeII fitting with a pre-set bindingSite, tested against
%%the SchemeII model with a Kdiss = 300, fixed.
elseif (inputPrompt == 4)
    indices = find(DNALength > bindingSite);    
    expr = ['vMaxATP*DNA/((A_uni/(DNALength-' num2str(bindingSite) ')) + K_diss + DNA)'];
    IIft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','A_uni','K_diss'});
    [IIfitModel,IIfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),IIft,'Start',[1 1 1],'Lower',[0 0 0]);
    confidenceInter = confint(IIfitModel);
    
    expr = ['vMaxATP*DNA/((A_uni/(DNALength-' num2str(bindingSite) ')) + 300 + DNA)'];
    restrictedIIft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','A_uni'});
    [restrictedIIfitModel,restrictedIIfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),restrictedIIft,'Start',[1 1],'Lower',[0 0]);
    confidenceInterRestricted = confint(restrictedIIfitModel);
    
    %Plotting of the fit to SchemeII from all relevant data.
    figure
    hold on
    subplot(2,1,1)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IIfitModel);
    titleString = ['Fitting data to SchemeII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(restrictedIIfitModel);
    titleString = ['Fitting data to restricted SchemeII with ' num2str(bindingSite) 'bp binding site and 300然 K_{diss}'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Plotting of the residual histograms.
    figure
    hold on
    subplot(2,1,1)
    plot(IIfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to SchemeII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    plot(restrictedIIfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to restricted SchemeII with ' num2str(bindingSite) 'bp binding site and 300然 K_{diss}'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Calculating the F-statistic of the model compared to the Null model
    %using the standard definition of the F-statistic for a restricted null
    %model. In the statistic below the Null model is always given by 1 and
    %the alternate hypothesis being the model tested.
    %FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2)
    %diff = vATP(:) - IIfitModel(DNALength(:),DNA(:));
    %SS2 = sum(diff(:).^2);
    SS2 = IIfitStatistics.sse;
    df2 = length(unique([DNALength,DNA], 'rows')) - 3;
    %diff = vATP(:) - restrictedIIfitModel(DNALength(:),DNA(:));
    %SS1 = sum(diff(:).^2);
    SS1 = restrictedIIfitStatistics.sse;
    df1 = length(unique([DNALength,DNA], 'rows')) - 2;
    FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2);
    pValue = 1-fcdf(FStatistic,(df1-df2),df2);
    
    %I want to display the final fitted parameters of the data and fits I have.
    disp('The models are fit to the SchemeII with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeII:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(IIfitModel.vMaxATP) '(' num2str(confidenceInter(1,1)) ',' num2str(confidenceInter(2,1)) ');'];
    disp(modelString);
    modelString = ['A_{uni} = ' num2str(IIfitModel.A_uni) '(' num2str(confidenceInter(1,2)) ',' num2str(confidenceInter(2,2)) ');'];
    disp(modelString);
    modelString = [' K_{diss} = ' num2str(IIfitModel.K_diss) '(' num2str(confidenceInter(1,3)) ',' num2str(confidenceInter(2,3)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IIfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('The models are fit to the Null Model SchemeI with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeI:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(restrictedIIfitModel.vMaxATP) '(' num2str(confidenceInterRestricted(1,1)) ',' num2str(confidenceInterRestricted(2,1)) ');'];
    disp(modelString);
    modelString = [' A_{uni} = ' num2str(restrictedIIfitModel.A_uni) '(' num2str(confidenceInterRestricted(1,2)) ',' num2str(confidenceInterRestricted(2,2)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(restrictedIIfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('By performing an F-test on the models and comparing them, the following is the p-value, i.e. the probability that null model is not rejected.');
    modelString = ['F-Statistic from F-test = ', num2str(FStatistic)];
    disp(modelString);
    modelString = ['p-Value from F-test = ', num2str(pValue)];
    disp(modelString);
    
%%This is the AlternateScheme(V) fitting with a [DNA] concentration that is calculated
%%as the effective binding site conc proportional to [DNA * DNALength].
%%vATP = (vMaxATP' * (DNA * DNALength))/(K'_ATP + (DNA * DNALength))
elseif (inputPrompt == 5)
    indices = find(DNALength > bindingSite);    
    expr = ['vMaxATP*(DNA*(DNALength-' num2str(bindingSite) '))/(Keff_ATP + DNA * (DNALength-' num2str(bindingSite) '))'];
    Vft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','Keff_ATP'});
    [VfitModel,VfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),Vft,'Start',[1 1],'Lower',[0 0]);
    confidenceInter = confint(VfitModel);
    
    expr = ['vMaxLength*(DNALength-' num2str(bindingSite) ')*DNA/((K_Length + (DNALength-' num2str(bindingSite) '))*(K_ATP + DNA))'];
    IIIft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxLength','K_Length','K_ATP'});
    [IIIfitModel,IIIfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),IIIft,'Start',[1 1 1],'Lower',[0 0 0]);
    confidenceInterIII = confint(IIIfitModel);
    
    %Plotting of the fit to SchemeV from all relevant data.
    figure
    hold on
    subplot(2,1,1)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(VfitModel);
    titleString = ['Fitting data to SchemeV with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IIIfitModel);
    titleString = ['Fitting data to SchemeIII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Plotting of the residual histograms.
    figure
    hold on
    subplot(2,1,1)
    plot(VfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to SchemeV with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    plot(IIIfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to SchemeIII with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Calculating the F-statistic of the model compared to the Null model
    %using the standard definition of the F-statistic for a restricted null
    %model. In the statistic below the Null model is always given by 1 and
    %the alternate hypothesis being the model tested.
    %FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2)
    %diff = vATP(:) - IIfitModel(DNALength(:),DNA(:));
    %SS2 = sum(diff(:).^2);
    SS2 = VfitStatistics.sse;
    df2 = length(unique([DNALength,DNA], 'rows')) - 2;
    %diff = vATP(:) - restrictedIIfitModel(DNALength(:),DNA(:));
    %SS1 = sum(diff(:).^2);
    SS1 = IIIfitStatistics.sse;
    df1 = length(unique([DNALength,DNA], 'rows')) - 3;
    %FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2);
    %pValue = 1-fcdf(FStatistic,(df1-df2),df2);
    
    %I want to display the final fitted parameters of the data and fits I have.
    disp('The models are fit to the SchemeV with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeV:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(VfitModel.vMaxATP) '(' num2str(confidenceInter(1,1)) ',' num2str(confidenceInter(2,1)) ');'];
    disp(modelString);
    modelString = [' Keff_{ATP} = ' num2str(VfitModel.Keff_ATP) '(' num2str(confidenceInter(1,2)) ',' num2str(confidenceInter(2,2)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(VfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('The models are fit to the SchemeIII with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeIII:';
    disp(modelString);
    modelString = ['v_{Max,Length} = ' num2str(IIIfitModel.vMaxLength) '(' num2str(confidenceInterIII(1,1)) ',' num2str(confidenceInterIII(2,1)) ');'];
    disp(modelString);
    modelString = ['K_{Length} = ' num2str(IIIfitModel.K_Length) '(' num2str(confidenceInterIII(1,2)) ',' num2str(confidenceInterIII(2,2)) ');'];
    disp(modelString);
    modelString = ['K_{ATP} = ' num2str(IIIfitModel.K_ATP) '(' num2str(confidenceInterIII(1,3)) ',' num2str(confidenceInterIII(2,3)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IIIfitStatistics.rsquare)];
    disp(modelString);
    
    %disp('By performing an F-test on the models and comparing them, the following is the p-value, i.e. the probability that null model is not rejected.');
    %modelString = ['F-Statistic from F-test = ', num2str(FStatistic)];
    %disp(modelString);
    %modelString = ['p-Value from F-test = ', num2str(pValue)];
    %disp(modelString);

%%This is the SchemeVI fitting with a pre-set bindingSite, tested against
%%the SchemeI Null Model.
elseif (inputPrompt == 6)
    indices = find(DNALength > bindingSite);    
    expr = ['vMaxLength*(DNALength-' num2str(bindingSite) ')*DNA/((K_Length + (DNALength-' num2str(bindingSite) '))*((A_uni/(DNALength-' num2str(bindingSite) ')) + K_diss + DNA))'];
    VIft = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxLength','K_Length','A_uni','K_diss'});
    [VIfitModel,VIfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),VIft,'Start',[1 1 1 1],'Lower',[0 0 0 0]);
    confidenceInter = confint(VIfitModel);
        
    expr = 'vMaxATP*DNA/(K_ATP + DNA) + 0*DNALength';
    Ift = fittype(expr,'dependent',{'vATP'},'independent',{'DNALength','DNA'},'coefficients',{'vMaxATP','K_ATP'});
    [IfitModel,IfitStatistics] = fit([DNALength(indices),DNA(indices)],vATP(indices),Ift,'Start',[1 1],'Lower',[0 0]);
    confidenceInterNull = confint(IfitModel);
    
    %Plotting of the fit to SchemeVI from all the data.
    figure
    hold on
    subplot(2,1,1)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(VIfitModel);
    titleString = ['Fitting data to SchemeVI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    scatter3(DNALength(indices),DNA(indices),vATP(indices));
    hold on; plot(IfitModel);
    titleString = ['Fitting data to Null Model SchemeI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Plotting of the residual histograms.
    figure
    hold on
    subplot(2,1,1)
    plot(VIfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to SchemeVI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    
    subplot(2,1,2)
    plot(IfitModel,[DNALength(indices) DNA(indices)], vATP(indices), 'Style', 'Residuals')
    titleString = ['Residuals of fit to Null Model SchemeI with ' num2str(bindingSite) 'bp binding site'];
    title(titleString);
    xlabel('DNA length (bp)');
    ylabel('[DNA] (然)');
    zlabel('vATP');
    hold off
    
    %Calculating the F-statistic of the model compared to the Null model
    %using the standard definition of the F-statistic for a restricted null
    %model. In the statistic below the Null model is always given by 1 and
    %the alternate hypothesis being the model tested.
    %FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2)
    %diff = vATP(:) - IIfitModel(DNALength(:),DNA(:));
    %SS2 = sum(diff(:).^2);
    SS2 = VIfitStatistics.sse;
    df2 = length(unique([DNALength,DNA], 'rows')) - 4;
    %diff = vATP(:) - IfitModel(DNALength(:),DNA(:));
    %SS1 = sum(diff(:).^2);
    SS1 = IfitStatistics.sse;
    df1 = length(unique([DNALength,DNA], 'rows')) - 2;
    FStatistic = ((SS1 - SS2)/(df1 - df2))/(SS2/df2);
    pValue = 1-fcdf(FStatistic,(df1-df2),df2);
    
    %I want to display the final fitted parameters of the data and fits I have.
    disp('The models are fit to the SchemeVI with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeVI:';
    disp(modelString);
    modelString = ['v_{Max,Length} = ' num2str(VIfitModel.vMaxLength) '(' num2str(confidenceInter(1,1)) ',' num2str(confidenceInter(2,1)) ');'];
    disp(modelString);
    modelString = ['K_{Length} = ' num2str(VIfitModel.K_Length) '(' num2str(confidenceInter(1,2)) ',' num2str(confidenceInter(2,2)) ');'];
    disp(modelString);
    modelString = ['A_{uni} = ' num2str(VIfitModel.A_uni) '(' num2str(confidenceInter(1,3)) ',' num2str(confidenceInter(2,3)) ');'];
    disp(modelString);
    modelString = ['K_{diss} = ' num2str(VIfitModel.K_diss) '(' num2str(confidenceInter(1,4)) ',' num2str(confidenceInter(2,4)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(VIfitStatistics.rsquare)];
    disp(modelString);
    
    disp('The models are fit to the Null Model SchemeI with the following parameters. The 95% confidence intervals and the R^{2} of the fit are also printed.');
    modelString = 'SchemeI:';
    disp(modelString);
    modelString = ['v_{Max,ATP} = ' num2str(IfitModel.vMaxATP) '(' num2str(confidenceInterNull(1,1)) ',' num2str(confidenceInterNull(2,1)) ');'];
    disp(modelString);
    modelString = [' K_{ATP} = ' num2str(IfitModel.K_ATP) '(' num2str(confidenceInterNull(1,2)) ',' num2str(confidenceInterNull(2,2)) ');'];
    disp(modelString);
    modelString = ['R^{2} = ', num2str(IfitStatistics.adjrsquare)];
    disp(modelString);
    
    disp('By performing an F-test on the models and comparing them, the following is the p-value, i.e. the probability that null model is not rejected.');
    modelString = ['F-Statistic from F-test = ', num2str(FStatistic)];
    disp(modelString);
    modelString = ['p-Value from F-test = ', num2str(pValue)];
    disp(modelString);
      
%%%Input check.
elseif (inputPrompt ~= (2 || 3 || 4 || 5 || 6))
    disp('You have to pick either 2 or 3 or 4 or 5 or 6!');
end
