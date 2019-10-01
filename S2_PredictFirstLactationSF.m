%% 1ST LACTATION vs RESILIENCE RANKING --> ONLY LACTATION SF
clear variables
close all
clc

% load data
load('D1_rankeddata.mat'); load('D3_crossvalidationDS.mat')
clear cd OUT2_1 OUT2_2 OUT2_3 OUT2_4 OUT2_5 SUM CL CR CLSF3 DAY_LT CLSF1 CLSF_LT cows_LT
Robot = {'DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY'};
Robot = double(contains(Robot,'DLV'));
Origin = {'KUL','KUL','KUL','KUL','KUL','KUL','KUL','KUL','KUL','KUL','RAFT','KUL','KUL','KUL','KUL','KUL','KUL','KUL','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT'};
Origin = double(contains(Origin,'KUL'));


% select SF and correct for outliers
SF = [17 24:52];
ind = find(CLSF2.Lac == 1 & CLSF2.NMilk > 200);
X = CLSF2(ind,[1 3 SF]);
ind = ismember(CRall{:,[1 3]},X{:,[1 2]},'rows');
Y = CRall(ind,[1 3 8 5]);
for iii = 3:length(X{1,:})
    X{find(isnan(X{:,iii})==1),iii}=0;  % replace NaN with average value = 0
end
clear iii i ii ind 


% PCA- check for underlying latent structure - not present
[~,~,PCAlatent,~,PCAexplained] = pca(X{:,:});
PCAexplained(:,2) = cumsum(PCAexplained(:,1));
clear PCAlatent PCAexplained iii


% Pearson linear correlations + figure
rho = array2table((1:length(SF))','VariableNames',{'IDtest'});
rho.SFName(:,1) = CLSF2.Properties.VariableNames(SF)';
for i = 1:27
    x = X(X.FarmID == i,3:end);
    y = sortrows(Y(Y.FarmID == i,:),2); y = (y{:,3}-1)./(length(y{:,1})-1);
        
    for ii = 1:length(x{1,:})
        rho{ii,i+2} = corr(x{:,ii},y,'Rows','complete');
    end
    
    F = sprintf('Farm_%d',i);
    rho.Properties.VariableNames{i+2} = F;
end
% add averages
rho.AvCor(:,1) = mean(rho{:,3:i+2},2);
rho.StdCor(:,1) = std(rho{:,3:i+2},[],2);
rho.LB(:,1) = rho.AvCor(:,1)-1.96*rho.StdCor(:,1);
rho.UB(:,1) = rho.AvCor(:,1)+1.96*rho.StdCor(:,1);
rho.AvUK(:,1) = mean(rho{:,2+find(Origin==0)},2);
rho.AvBE(:,1) = mean(rho{:,2+find(Origin==1)},2);
% plot results
figure; hold on; box on; title('Pearsons linear correlations first lactation feature vs resilience RANK per farm, blue = UK, purple = BE'); LW = 'LineWidth'; 
xlabel('Sensor feature'); ylabel('Correlation coefficient rho'); xlim([0.8 30.2 ]); ylim([-0.55 0.55])
xticks(1:length(SF)); xticklabels(CLSF2.Properties.VariableNames(SF)); xtickangle(45)
h1 = area(1:30,rho.UB(:,1),-0.4); h1(1).FaceColor = [0.8 0.8 0.8]; h1(1).EdgeColor = [1 1 1];
h2 = area(1:30,rho.LB(:,1),-0.4); h2(1).FaceColor = [1 1 1 ]; h2(1).EdgeColor = [1 1 1];
plot(1:30,zeros(1,30),'r--','LineWidth',2.5)
r = 0.1:0.03:0.9;
% FARM
r3 = 0.1:0.082:0.9; % for UK, n = 10
r4 = 0.1:0.05:0.9; % for BE, n = 17
COL3 = [r3'.^6.8 r3'.^1.7 r3']; T1 = 1;
COL4 = [r4' r4'.^6.5 r4'.^1.2]; T2 = 1;
for i = 1:27
    F = sprintf('Farm_%d',i);    
    if Origin(i) == 0 % Robot(i) == 1
        plot(rho.(F),'o-',LW, 1,'Color',COL3(T1,:),'MarkerFaceColor',COL3(T1,:),'MarkerSize',5);
        T1 = T1+1;
    else
        plot(rho.(F),'o-',LW, 1,'Color',COL4(T2,:),'MarkerFaceColor',COL4(T2,:),'MarkerSize',5);
        T2 = T2+1;
    end
end
plot(rho.AvCor,'k','LineWidth',3.5)
plot(rho.AvUK,'Color',[0 0.75 0.85],LW,3.5) % blueish
plot(rho.AvBE,'Color',[0.75, 0, 0.75],LW,3.5) % redish
rho{31,3:end} = max(rho{1:length(SF)-1,3:end}); rho{31,2} = {'max'};
rho{32,3:end} = min(rho{1:length(SF)-1,3:end}); rho{32,2} = {'min'};
clear ans COL3 COL4 F h1 h2 i ii LW r r3 r4 T1 T2 x y


% multilinear regression
b = array2table((1:length(SF))','VariableNames',{'IDtest'});
b.SFName(:,1) = X.Properties.VariableNames(3:end)';
b.IDtest(end+1:end+4) = length(SF)+1:length(SF)+4'; 
b.SFName(length(SF)+1)= {'Intercept'}; b.SFName(length(SF)+2) = {'R2'};b.SFName(length(SF)+3) = {'RMSE'};
b.SFName(length(SF)+4) = {'Nterms'};

for i = 1:27
    
    x = X(find(X.FarmID == i),3:end);
    
    for iii = 1:length(X{1,:})
        X{find(isnan(X{:,iii})==1),iii}=0;  % replace NaN with average value = 0
    end
    
    y = sortrows(Y(Y.FarmID == i,:),2); y = (y{:,3}-1)./(max(y{:,3})-1); % y and scale y
    
    [STPWSb,~,~,STPWSinmodel,STPWSstats] = stepwisefit(x{:,:},y, 'penter',0.2,'premove',0.2,'maxiter',200,'display','off');
    
    b{STPWSinmodel',i+2} = STPWSb(STPWSinmodel');
    b{length(SF)+1,i+2} = STPWSstats.intercept;
    R = 1-(STPWSstats.SSresid./STPWSstats.SStotal);
    n = length(y); k = sum(STPWSinmodel);
    b{length(SF)+2,i+2} = 1-(((1-R)*(n-1))/(n-k-1));
    b{length(SF)+3,i+2} = STPWSstats.rmse;
    b{length(SF)+4,i+2} = sum(STPWSinmodel);
    
    Yfitted = b{1:length(SF)+1,i+2}'*[x{:,:} ones(length(x{:,1}),1)]';
    
    figure(12); subplot(4,7,i); hold on; ylabel('RR_{pred}');xlabel('RR'); box on; xlim([0 1]); ylim([0 1])
    plot(y, Yfitted, '.','MarkerSize',10);
    plot([0 1],[0 1],'r','LineWidth',2)
    title(['F',num2str(i),', Rs=' num2str(round(b{length(SF)+2,i+2},2)) ', n=' num2str(length(find(b{1:length(SF)+1,i+2})))])
        
    F = sprintf('Farm_%d',i);
    b.Properties.VariableNames{i+2} = F;
end
clear x y Yfitted STPWSb STPWSstats F i ii iii k n R


% cross validation of the linear regression
GOF = []; 
for ii = 1:10              % cross validation 1 to 10
    CV = sprintf('CV_%d',ii);
    
    GOF.(CV) = array2table((1:25)','VariableNames',{'IDtest'});
    GOF.(CV).Name(:,1) = {'R2fit','R2pred','RMSEC','RMSECV','RMSEvalCor','BB','BM','BW','MB','MM','MW','WB','WM','WW','BBc','BMc','BWc','MBc','MMc','MWc','WBc','WMc','WWc','FirstCor','HighCor'}';
    T = 1;

    for i = 1:27
        
        SF2 = SF(find(b{1:length(SF),i+2}));   % select SF for specific farm
        
        x = X(find(X.FarmID == i),2+find(b{1:length(SF),i+2}));

        testindx = find(ismember(find(X.FarmID == i),test(:,ii))==1);   % select indices of testset
        traiindx = find(ismember(find(X.FarmID == i),training(:,ii))==1);   % select indices of training set

        y = sortrows(Y(Y.FarmID == i,:),2);     % select y
        lastlac = y;                            % store info for last = first selection
        y = (y{:,3}-1)./(max(y{:,3})-1);        % scale y
        
        Yreg = y(traiindx);         % ytraining
        Xreg = x{traiindx,:};       % xtraining
        
        VN = [CLSF2.Properties.VariableNames(SF2)' ; 'RRank'];
        F = sprintf('Farm_%d',i);

        ModVal = fitlm(Xreg,Yreg,'VarNames',VN);  % fit model
        
        TrainFit = ModVal.Fitted';
        A = min(TrainFit); B = max(TrainFit); 

        Yval = y(testindx);         % yvalidation
        Xval = x{testindx,:};       % xvalidation
        lastlacval = lastlac(testindx,:); % last lactation of the testset
        
        fitted = (ModVal.Coefficients.Estimate' * [ones(length(Xval(:,1)),1) Xval]')';
        resid = Yval - fitted;
        RMSEpred = sqrt(nanmean(resid.^2));
        RMSEfit = sqrt(nanmean(ModVal.Residuals.Raw.^2));
        SSE = sum(resid.^2);
        SSTO = sum((Yval-mean(Yval)).^2);
                
        fittedcorrected = (fitted-A)./(B-A);        % on validation dataset
        fittedcorrected(:,2) = Yval;
        residcorrected = Yval - fittedcorrected(:,1);        % on training dataset
        RMSEfitcorrected = sqrt(nanmean(residcorrected.^2));       % 
       
        lastlacval.Yval(:,1) = fittedcorrected(:,2);     % the actual RR (corrected)
        lastlacval.Ypred(:,1) = fittedcorrected(:,1);     % the actual RR (corrected)
        lastlacval.Corr(:,1) = 0;                       % to fill in with 1 if ocrrectly predicted
        lastlacval.Corr((lastlacval.Yval <= 0.33 & lastlacval.Ypred <= 0.33)|((lastlacval.Yval > 0.33 & lastlacval.Yval < 0.66) & (lastlacval.Ypred > 0.33 & lastlacval.Ypred < 0.66))|(lastlacval.Yval >= 0.66 & lastlacval.Ypred >= 0.66)) = 1;
        
        % prediction accuracy
        BB = length(find(Yval <= 0.33 & fitted <= 0.33));                       % true = best, pred = best
        BM = length(find(Yval <= 0.33 & (fitted > 0.33 & fitted < 0.66)));      % true = best, pred = middle
        BW = length(find(Yval <= 0.33 & fitted >= 0.66));                       % true = best, pred = worst
        MB = length(find((Yval > 0.33 & Yval < 0.66) & fitted <= 0.33));        % true = middle, pred = best
        MM = length(find((Yval > 0.33 & Yval < 0.66) & (fitted > 0.33 & fitted < 0.66)));      % true = middle, pred = middle
        MW = length(find((Yval > 0.33 & Yval < 0.66) & fitted >= 0.66));                       % true = middle, pred = worst
        WB = length(find(Yval > 0.66 & fitted <= 0.33));                       % true = worst, pred = best
        WM = length(find(Yval > 0.66 & (fitted > 0.33 & fitted < 0.66)));      % true = worst, pred = middle
        WW = length(find(Yval > 0.66 & fitted >= 0.66));                       % true = worst, pred = worst
        
        % prediction accuracy after correction
        BBc = length(find(Yval <= 0.33 & fittedcorrected(:,1) <= 0.33));                       % true = best, pred = best
        BMc = length(find(Yval <= 0.33 & (fittedcorrected(:,1) > 0.33 & fittedcorrected(:,1) < 0.66)));      % true = best, pred = middle
        BWc = length(find(Yval <= 0.33 & fittedcorrected(:,1) >= 0.66));                       % true = best, pred = worst
        MBc = length(find((Yval > 0.33 & Yval < 0.66) & fittedcorrected(:,1) <= 0.33));        % true = middle, pred = best
        MMc = length(find((Yval > 0.33 & Yval < 0.66) & (fittedcorrected(:,1) > 0.33 & fittedcorrected(:,1) < 0.66)));      % true = middle, pred = middle
        MWc = length(find((Yval > 0.33 & Yval < 0.66) & fittedcorrected(:,1) >= 0.66));                       % true = middle, pred = worst
        WBc = length(find(Yval > 0.66 & fittedcorrected(:,1) <= 0.33));                       % true = worst, pred = best
        WMc = length(find(Yval > 0.66 & (fittedcorrected(:,1) > 0.33 & fittedcorrected(:,1) < 0.66)));      % true = worst, pred = middle
        WWc = length(find(Yval > 0.66 & fittedcorrected(:,1) >= 0.66));                       % true = worst, pred = worst

                
        % store GOF measures 
        GOF.(CV).(F)(1,1) = ModVal.Rsquared.Ordinary(1,1);        % fitted R2
        GOF.(CV).(F)(2,1) = 1-SSE/SSTO;        % fitted R2
        GOF.(CV).(F)(3,1) = RMSEfit;             % RMSE fit
        GOF.(CV).(F)(4,1) = RMSEpred;
        GOF.(CV).(F)(5,1) = RMSEfitcorrected;       % RMSEfitcorr
        GOF.(CV).(F)(6,1) = BB;
        GOF.(CV).(F)(7,1) = BM;
        GOF.(CV).(F)(8,1) = BW;
        GOF.(CV).(F)(9,1) = MB;
        GOF.(CV).(F)(10,1) = MM;
        GOF.(CV).(F)(11,1) = MW;
        GOF.(CV).(F)(12,1) = WB;
        GOF.(CV).(F)(13,1) = WM;
        GOF.(CV).(F)(14,1) = WW;
        GOF.(CV).(F)(15,1) = BBc;
        GOF.(CV).(F)(16,1) = BMc;
        GOF.(CV).(F)(17,1) = BWc;
        GOF.(CV).(F)(18,1) = MBc;
        GOF.(CV).(F)(19,1) = MMc;
        GOF.(CV).(F)(20,1) = MWc;
        GOF.(CV).(F)(21,1) = WBc;
        GOF.(CV).(F)(22,1) = WMc;
        GOF.(CV).(F)(23,1) = WWc;
        GOF.(CV).(F)(24,1) = sum(lastlacval.Corr == 1 & lastlacval.LastLac == 1);
        GOF.(CV).(F)(25,1) = sum(lastlacval.Corr == 1 & lastlacval.LastLac > 1);


        T = T+1;
    end
    clear SSE SSTO RMSEfit RMSEpred PREDBEST PREDWORST WORST BEST iii Xreg Xval Yreg Yval F VN traiindx testindx SSE SSTO x y
end
clear T  resid i ii CV ans ALFIRST RMSE RMSEBEST ModVal SF2 fitted BB BM BW MM MB MW WB WM WW STPWSb STPWSinmodel STPWSstats Yfitted


% Overall performance over the CV
GOF.ALL = array2table((1:25)','VariableNames',{'IDtest'});
GOF.ALL.Name(:,1) = {'R2fit','R2pred','RMSEC','RMSECV','RMSECVcor','BB','BM','BW','MB','MM','MW','WB','WM','WW','BBc','BMc','BWc','MBc','MMc','MWc','WBc','WMc','WWc','FirstCor','HighCor'}';

for ii = 1:10
    CV = sprintf('CV_%d',ii);
    
    AV_R2fit(ii,:) = GOF.(CV)(1,3:29);
    AV_R2pred(ii,:) = GOF.(CV)(2,3:29);
    AV_RMSEfit(ii,:) = GOF.(CV)(3,3:29);
    AV_RMSEpred(ii,:) = GOF.(CV)(4,3:29);
    AV_RMSEpredcorr(ii,:) = GOF.(CV)(5,3:29);
    N_BB(ii,:) = GOF.(CV)(6,3:29);
    N_BM(ii,:) = GOF.(CV)(7,3:29);
    N_BW(ii,:) = GOF.(CV)(8,3:29);
    N_MB(ii,:) = GOF.(CV)(9,3:29);
    N_MM(ii,:) = GOF.(CV)(10,3:29);
    N_MW(ii,:) = GOF.(CV)(11,3:29);
    N_WB(ii,:) = GOF.(CV)(12,3:29);
    N_WM(ii,:) = GOF.(CV)(13,3:29);
    N_WW(ii,:) = GOF.(CV)(14,3:29);
    N_BBc(ii,:) = GOF.(CV)(15,3:29);
    N_BMc(ii,:) = GOF.(CV)(16,3:29);
    N_BWc(ii,:) = GOF.(CV)(17,3:29);
    N_MBc(ii,:) = GOF.(CV)(18,3:29);
    N_MMc(ii,:) = GOF.(CV)(19,3:29);
    N_MWc(ii,:) = GOF.(CV)(20,3:29);
    N_WBc(ii,:) = GOF.(CV)(21,3:29);
    N_WMc(ii,:) = GOF.(CV)(22,3:29);
    N_WWc(ii,:) = GOF.(CV)(23,3:29);
    N_FirstCor(ii,:) = GOF.(CV)(24,3:29);
    N_HighCor(ii,:) = GOF.(CV)(25,3:29);
end
for i = 1:27
    AV_R2pred{find(AV_R2pred{:,i}<0),i} = 0;
end 

GOF.ALL{1,3:29} = nanmean(AV_R2fit{:,:});
GOF.ALL{2,3:29} = nanmean(AV_R2pred{:,:}); 
GOF.ALL{3,3:29} = nanmean(AV_RMSEfit{:,:});
GOF.ALL{4,3:29} = nanmean(AV_RMSEpred{:,:});
GOF.ALL{5,3:29} = nanmean(AV_RMSEpredcorr{:,:});
GOF.ALL{6,3:29} = nansum(N_BB{:,:});
GOF.ALL{7,3:29} = nansum(N_BM{:,:});
GOF.ALL{8,3:29} = nansum(N_BW{:,:});
GOF.ALL{9,3:29} = nansum(N_MB{:,:});
GOF.ALL{10,3:29} = nansum(N_MM{:,:});
GOF.ALL{11,3:29} = nansum(N_MW{:,:});
GOF.ALL{12,3:29} = nansum(N_WB{:,:});
GOF.ALL{13,3:29} = nansum(N_WM{:,:});
GOF.ALL{14,3:29} = nansum(N_WW{:,:});
GOF.ALL{15,3:29} = nansum(N_BBc{:,:});
GOF.ALL{16,3:29} = nansum(N_BMc{:,:});
GOF.ALL{17,3:29} = nansum(N_BWc{:,:});
GOF.ALL{18,3:29} = nansum(N_MBc{:,:});
GOF.ALL{19,3:29} = nansum(N_MMc{:,:});
GOF.ALL{20,3:29} = nansum(N_MWc{:,:});
GOF.ALL{21,3:29} = nansum(N_WBc{:,:});
GOF.ALL{22,3:29} = nansum(N_WMc{:,:});
GOF.ALL{23,3:29} = nansum(N_WWc{:,:});

GOF.ALL{34,3:29} = nansum(N_FirstCor{:,:}); GOF.ALL{34,2} = {'FirstCor'};       
GOF.ALL{35,3:29} = nansum(N_HighCor{:,:});  GOF.ALL{35,2} = {'HighCor'};       


GOF.ALL{24,3:29} = sum(GOF.ALL{6:8,3:29}); GOF.ALL{24,2} = {'B'};       
GOF.ALL{25,3:29} = sum(GOF.ALL{9:11,3:29}); GOF.ALL{25,2} = {'M'}; 
GOF.ALL{26,3:29} = sum(GOF.ALL{12:14,3:29}); GOF.ALL{26,2} = {'W'};
GOF.ALL{27,3:29} = sum(GOF.ALL{15:17,3:29}); GOF.ALL{27,2} = {'Bc'};
GOF.ALL{28,3:29} = sum(GOF.ALL{18:20,3:29}); GOF.ALL{28,2} = {'Mc'}; 
GOF.ALL{29,3:29} = sum(GOF.ALL{21:23,3:29}); GOF.ALL{29,2} = {'Wc'};

% true positive and complete opposite %
GOF.ALL{30,3:29} = (GOF.ALL{5,3:29}+GOF.ALL{9,3:29}+GOF.ALL{13,3:29})./sum(GOF.ALL{24:26,3:29}); GOF.ALL{30,2} = {'TRUEPRED'};  GOF.ALL{30,1} = 30;  % 
GOF.ALL{31,3:29} = (GOF.ALL{15,3:29}+GOF.ALL{19,3:29}+GOF.ALL{23,3:29})./sum(GOF.ALL{27:29,3:29}); GOF.ALL{31,2} = {'TRUEPREDc'};  GOF.ALL{31,1} = 31;  % 
GOF.ALL{32,3:29} = (GOF.ALL{8,3:29}+GOF.ALL{12,3:29})./sum(GOF.ALL{24:26,3:29}); GOF.ALL{32,2} = {'OPPOSITE'};  GOF.ALL{32,1} = 32;
GOF.ALL{33,3:29} = (GOF.ALL{17,3:29}+GOF.ALL{21,3:29})./sum(GOF.ALL{24:26,3:29}); GOF.ALL{33,2} = {'OPPOSITEc'};  GOF.ALL{33,1} = 33;

% correctly predicted and in first lactation
GOF.ALL{36,3:29} = GOF.ALL{34,3:29}./(GOF.ALL{15,3:29}+GOF.ALL{19,3:29}+GOF.ALL{23,3:29});  % percentage correctly predicted in first lactation
GOF.ALL{37,3:29} = GOF.ALL{35,3:29}./(GOF.ALL{15,3:29}+GOF.ALL{19,3:29}+GOF.ALL{23,3:29});  % percentage correctly predicted in higher lactation

% total number of 1st lactation = last lactation cows in test sets
for ii = 1:10
    for i = 1:27
        testindx = find(ismember(find(X.FarmID == i),test(:,ii))==1);   % select indices of testset

        y = sortrows(Y(Y.FarmID == i,:),2);     % select y
        
        Nfirstlact = y(testindx,:);         % ytraining
        Testsetlact1(i,ii) = length(find(Nfirstlact.LastLac == 1));
        Testsetlact2(i,ii) = length(find(Nfirstlact.LastLac > 1));
    end
end
percfirsttestset = Testsetlact1./(Testsetlact1+Testsetlact2);
percfirsttestset(:,11) = mean(percfirsttestset(:,1:10),2)*100;

% per validation % first lactation and first+correctly predicted
perc_corr_first = (N_FirstCor{:,:}')./Testsetlact1;
perc_corr_first(:,11) = nanmean(perc_corr_first(:,1:10),2);

clear N_TH N_TH_PH N_TL N_TL_PH i  ii fitted CV
clear N_BB N_BM N_BW N_MM N_MW N_MB N_WW N_WB N_WM N_BBc N_BMc N_BWc N_MMc N_MWc N_MBc N_WWc N_WBc N_WMc
clear WBc WMc WWc MBc MMc MWc BBc BMc BWc B A TrainFit  RMSEfitcorr residcorrected fittedcorrected AV_R2fit AV_R2pred AV_RMSEpred AV_RMSEpredcorr AV_RMSEfit RMSEfitcorrected



%% 1ST LACTATION vs RESILIENCE RANKING --> LACTATION & ACTIVITY SF
clear variables
close all
clc

% load data
load('D1_rankeddata.mat'); load('D3_crossvalidationDS.mat'); load('D2_ActivityFeatures.mat')
CLSF2.Properties.VariableNames = [CLSF2.Properties.VariableNames(1:86) CLSF_ACT.Bentley.Properties.VariableNames(7:48)];
clear cd OUT2_1 OUT2_2 OUT2_3 OUT2_4 OUT2_5 SUM CL CR CLSF3 DAY_LT CLSF1 CLSF_LT cows_LT CLSF_ACT FN
Robot = {'DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','DLV','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY','LELY'};
Robot = double(contains(Robot,'DLV'));
Origin = {'KUL','KUL','KUL','KUL','KUL','KUL','KUL','KUL','KUL','KUL','RAFT','KUL','KUL','KUL','KUL','KUL','KUL','KUL','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT','RAFT'};
Origin = double(contains(Origin,'KUL'));
Farms = [12 13 16 17 18 19 21 22 23 24 25 26 27];         % all farm ID of farms with activity data = 13 farms
Robot = Robot(Farms); Origin = Origin(Farms);

% select SF and correct for outliers
SF = [17 24:52 89 92 93 95 96 97 103 104 105 106 109 113 114 115 118 ]; % 46 features

ind = find(CLSF2.Lac == 1 & CLSF2.NMilk > 200 & ismember(CLSF2.FarmID,Farms)==1);
X = CLSF2(ind,[1 3 SF]);
ind = ismember(CRall{:,[1 3]},X{:,[1 2]},'rows');
Y = CRall(ind,[1 3 8 5]);
for iii = 3:length(X{1,:})
    X{find(isnan(X{:,iii})==1),iii}=0;  % replace NaN with average value = 0
end
clear iii i ii ind 


% PCA- check for underlying latent structure - not present
[~,~,PCAlatent,~,PCAexplained] = pca(X{:,3:end});
PCAexplained(:,2) = cumsum(PCAexplained(:,1));
clear PCAlatent PCAexplained iii


% Pearson linear correlations + figure
rho = array2table((1:length(SF))','VariableNames',{'IDtest'});
rho.SFName(:,1) = CLSF2.Properties.VariableNames(SF)';
T = 1;
for i = Farms
    x = X(X.FarmID == i,3:end);
    y = sortrows(Y(Y.FarmID == i,:),2); y = (y{:,3}-1)./(length(y{:,1})-1);
        
    for ii = 1:length(x{1,:})
        rho{ii,T+2} = corr(x{:,ii},y,'Rows','complete');
    end
    
    F = sprintf('Farm_%d',Farms(T));
    rho.Properties.VariableNames{T+2} = F;
    T = T+1;
end
% add averages
rho.AvCor(:,1) = mean(rho{:,3:T+1},2);
rho.StdCor(:,1) = std(rho{:,3:T+1},[],2);
rho.LB(:,1) = rho.AvCor(:,1)-1.96*rho.StdCor(:,1);
rho.UB(:,1) = rho.AvCor(:,1)+1.96*rho.StdCor(:,1);
rho.AvUK(:,1) = mean(rho{:,2+find(Origin==0)},2);
rho.AvBE(:,1) = mean(rho{:,2+find(Origin==1)},2);
% plot results
figure; hold on; box on; title('Pearsons linear correlations first lactation feature vs resilience RANK per farm, blue = UK, purple = BE'); LW = 'LineWidth'; 
xlabel('Sensor feature'); ylabel('Correlation coefficient rho'); xlim([0.8 45.2 ]); ylim([-0.55 0.55])
xticks(1:length(SF)); xticklabels(CLSF2.Properties.VariableNames(SF)); xtickangle(45)
h1 = area(1:45,rho.UB(:,1),-0.4); h1(1).FaceColor = [0.8 0.8 0.8]; h1(1).EdgeColor = [1 1 1];
h2 = area(1:45,rho.LB(:,1),-0.4); h2(1).FaceColor = [1 1 1 ]; h2(1).EdgeColor = [1 1 1];
plot(1:45,zeros(1,45),'r--','LineWidth',2.5)
% FARM
r3 = 0.1:0.11:0.9; % for BE, n = 8
r4 = 0.1:0.17:0.9; % for UK, n = 5
COL3 = [r3'.^6.8 r3'.^1.7 r3']; T1 = 1;  % BE
COL4 = [r4' r4'.^6.5 r4'.^1.2]; T2 = 1;  % UK
T = 1;
for i = Farms
    F = sprintf('Farm_%d',i);    
    if Origin(T) == 0   % UK farms
        plot(rho.(F),'o-',LW, 1,'Color',COL3(T1,:),'MarkerFaceColor',COL3(T1,:),'MarkerSize',5);
        T1 = T1+1;
    else  % orgin is BE
        plot(rho.(F),'o-',LW, 1,'Color',COL4(T2,:),'MarkerFaceColor',COL4(T2,:),'MarkerSize',5);
        T2 = T2+1;
    end
    T = T+1;
end
plot(rho.AvCor,'k','LineWidth',3.5)
plot(rho.AvUK,'Color',[0 0.75 0.85],LW,3.5) % blueish
plot(rho.AvBE,'Color',[0.75, 0, 0.75],LW,3.5) % redish
rho{46,3:end} = max(rho{1:length(SF)-1,3:end}); rho{46,2} = {'max'};
rho{47,3:end} = min(rho{1:length(SF)-1,3:end}); rho{47,2} = {'min'};
clear ans COL3 COL4 F h1 h2 i ii LW r r3 r4 T1 T2 x y

% plot only activity features
figure; hold on; box on; title('Pearsons linear correlations first lactation activity features vs resilience RANK per farm, blue = UK, purple = BE'); LW = 'LineWidth'; 
xlabel('Sensor feature'); ylabel('Correlation coefficient rho'); xlim([0.8 15.2 ]); ylim([-0.55 0.55])
xticks(1:length(SF)); xticklabels(CLSF2.Properties.VariableNames(SF(31:45))); xtickangle(45)
h1 = area(1:15,rho.UB(31:45,1),-0.4); h1(1).FaceColor = [0.8 0.8 0.8]; h1(1).EdgeColor = [1 1 1];
h2 = area(1:15,rho.LB(31:45,1),-0.4); h2(1).FaceColor = [1 1 1 ]; h2(1).EdgeColor = [1 1 1];
plot(1:15,zeros(1,15),'r--','LineWidth',2.5)
% FARM
r3 = 0.1:0.11:0.9; % for BE, n = 8
r4 = 0.1:0.17:0.9; % for UK, n = 5
COL3 = [r3'.^6.8 r3'.^1.7 r3']; T1 = 1;  % BE
COL4 = [r4' r4'.^6.5 r4'.^1.2]; T2 = 1;  % UK
T = 1;
for i = Farms
    F = sprintf('Farm_%d',i);    
    if Origin(T) == 0   % UK farms
        plot(rho.(F)(31:45),'o-',LW, 1,'Color',COL3(T1,:),'MarkerFaceColor',COL3(T1,:),'MarkerSize',5);
        T1 = T1+1;
    else  % orgin is BE
        plot(rho.(F)(31:45),'o-',LW, 1,'Color',COL4(T2,:),'MarkerFaceColor',COL4(T2,:),'MarkerSize',5);
        T2 = T2+1;
    end
    T = T+1;
end
plot(rho.AvCor(31:45),'k','LineWidth',3.5)
plot(rho.AvUK(31:45),'Color',[0 0.75 0.85],LW,3.5) % blueish
plot(rho.AvBE(31:45),'Color',[0.75, 0, 0.75],LW,3.5) % redish
rho{46,3:end} = max(rho{1:length(SF)-1,3:end}); rho{46,2} = {'max'};
rho{47,3:end} = min(rho{1:length(SF)-1,3:end}); rho{47,2} = {'min'};
clear ans COL3 COL4 F h1 h2 i ii LW r r3 r4 T1 T2 x y




% multilinear regression
b = array2table((1:length(SF))','VariableNames',{'IDtest'});
b.SFName(:,1) = X.Properties.VariableNames(3:end)';
b.IDtest(end+1:end+4) = length(SF)+1:length(SF)+4'; 
b.SFName(length(SF)+1)= {'Intercept'}; b.SFName(length(SF)+2) = {'R2'};b.SFName(length(SF)+3) = {'RMSE'};
b.SFName(length(SF)+4) = {'Nterms'};
T = 1;
for i = Farms
    
    x = X(find(X.FarmID == i),3:end);
    
    for iii = 1:length(X{1,:})
        X{find(isnan(X{:,iii})==1),iii}=0;  % replace NaN with average value = 0
    end
    
    y = sortrows(Y(Y.FarmID == i,:),2); y = (y{:,3}-1)./(max(y{:,3})-1); % y and scale y
    
    [STPWSb,~,~,STPWSinmodel,STPWSstats] = stepwisefit(x{:,:},y, 'penter',0.2,'premove',0.2,'maxiter',200,'display','off');
    
    b{STPWSinmodel',T+2} = STPWSb(STPWSinmodel');
    b{length(SF)+1,T+2} = STPWSstats.intercept;
    R = 1-(STPWSstats.SSresid./STPWSstats.SStotal);
    n = length(y); k = sum(STPWSinmodel);
    b{length(SF)+2,T+2} = 1-(((1-R)*(n-1))/(n-k-1));
    b{length(SF)+3,T+2} = STPWSstats.rmse;
    b{length(SF)+4,T+2} = sum(STPWSinmodel);
    
    Yfitted = b{1:length(SF)+1,T+2}'*[x{:,:} ones(length(x{:,1}),1)]';
    
    figure(13); subplot(3,5,T); hold on; ylabel('RR_{pred}');xlabel('RR'); box on; xlim([0 1]); ylim([0 1])
    plot(y, Yfitted, '.','MarkerSize',10);
    plot([0 1],[0 1],'r','LineWidth',2)
    title(['F',num2str(i),', Rs=' num2str(round(b{length(SF)+2,T+2},2)) ', n=' num2str(length(find(b{1:length(SF)+1,T+2})))])
        
    F = sprintf('Farm_%d',i);
    b.Properties.VariableNames{T+2} = F;
    T = T+1;
end
clear x y Yfitted STPWSb STPWSstats F i ii iii k n R


% cross validation of the linear regression
GOF = []; 
for ii = 1:10              % cross validation 1 to 10
    CV = sprintf('CV_%d',ii);
    
    GOF.(CV) = array2table((1:25)','VariableNames',{'IDtest'});
    GOF.(CV).Name(:,1) = {'R2fit','R2pred','RMSEC','RMSECV','RMSEvalCor','BB','BM','BW','MB','MM','MW','WB','WM','WW','BBc','BMc','BWc','MBc','MMc','MWc','WBc','WMc','WWc','FirstCor','HighCor'}';
    T = 1;
    TT = 1;
    for i = Farms
        
        SF2 = SF(find(b{1:length(SF),TT+2}));   % select SF for specific farm
        
        x = X(find(X.FarmID == i),2+find(b{1:length(SF),TT+2}));

        testindx = find(ismember(find(X.FarmID == i),test(:,ii))==1);
        traiindx = find(ismember(find(X.FarmID == i),training(:,ii))==1);

        y = sortrows(Y(Y.FarmID == i,:),2);     % select y
        lastlac = y;                            % store info for last = first selection
        y = (y{:,3}-1)./(max(y{:,3})-1);        % scale y
        
        Yreg = y(traiindx);         % ytraining
        Xreg = x{traiindx,:};       % xtraining

        VN = [CLSF2.Properties.VariableNames(SF2)' ; 'RRank'];
        F = sprintf('Farm_%d',i);

        ModVal = fitlm(Xreg,Yreg,'VarNames',VN);  % fit model
        
        TrainFit = ModVal.Fitted';
        A = min(TrainFit); B = max(TrainFit); 

        Yval = y(testindx);         % yvalidation
        Xval = x{testindx,:};       % xvalidation
        lastlacval = lastlac(testindx,:); % last lactation of the cows in testset
        
        fitted = (ModVal.Coefficients.Estimate' * [ones(length(Xval(:,1)),1) Xval]')';
        resid = Yval - fitted;
        RMSEpred = sqrt(nanmean(resid.^2));
        RMSEfit = sqrt(nanmean(ModVal.Residuals.Raw.^2));
        SSE = sum(resid.^2);
        SSTO = sum((Yval-mean(Yval)).^2);
                
        fittedcorrected = (fitted-A)./(B-A);        % on validation dataset
        fittedcorrected(:,2) = Yval;
        residcorrected = Yval - fittedcorrected(:,1);        % on training dataset
        RMSEfitcorrected = sqrt(nanmean(residcorrected.^2));       % 
       
        lastlacval.Yval(:,1) = fittedcorrected(:,2);     % the actual RR (corrected)
        lastlacval.Ypred(:,1) = fittedcorrected(:,1);     % the actual RR (corrected)
        lastlacval.Corr(:,1) = 0;                       % to fill in with 1 if ocrrectly predicted
        lastlacval.Corr((lastlacval.Yval <= 0.33 & lastlacval.Ypred <= 0.33)|((lastlacval.Yval > 0.33 & lastlacval.Yval < 0.66) & (lastlacval.Ypred > 0.33 & lastlacval.Ypred < 0.66))|(lastlacval.Yval >= 0.66 & lastlacval.Ypred >= 0.66)) = 1;

   
        % prediction accuracy
        BB = length(find(Yval <= 0.33 & fitted <= 0.33));                       % true = best, pred = best
        BM = length(find(Yval <= 0.33 & (fitted > 0.33 & fitted < 0.66)));      % true = best, pred = middle
        BW = length(find(Yval <= 0.33 & fitted >= 0.75));                       % true = best, pred = worst
        MB = length(find((Yval > 0.33 & Yval < 0.66) & fitted <= 0.33));        % true = middle, pred = best
        MM = length(find((Yval > 0.33 & Yval < 0.66) & (fitted > 0.33 & fitted < 0.66)));      % true = middle, pred = middle
        MW = length(find((Yval > 0.33 & Yval < 0.66) & fitted >= 0.66));                       % true = middle, pred = worst
        WB = length(find(Yval > 0.66 & fitted <= 0.33));                       % true = worst, pred = best
        WM = length(find(Yval > 0.66 & (fitted > 0.33 & fitted < 0.66)));      % true = worst, pred = middle
        WW = length(find(Yval > 0.66 & fitted >= 0.66));                       % true = worst, pred = worst
        
        % prediction accuracy after correction
        BBc = length(find(Yval <= 0.33 & fittedcorrected(:,1) <= 0.33));                       % true = best, pred = best
        BMc = length(find(Yval <= 0.33 & (fittedcorrected(:,1) > 0.33 & fittedcorrected(:,1) < 0.66)));      % true = best, pred = middle
        BWc = length(find(Yval <= 0.33 & fittedcorrected(:,1) >= 0.75));                       % true = best, pred = worst
        MBc = length(find((Yval > 0.33 & Yval < 0.66) & fittedcorrected(:,1) <= 0.33));        % true = middle, pred = best
        MMc = length(find((Yval > 0.33 & Yval < 0.66) & (fittedcorrected(:,1) > 0.33 & fittedcorrected(:,1) < 0.66)));      % true = middle, pred = middle
        MWc = length(find((Yval > 0.33 & Yval < 0.66) & fittedcorrected(:,1) >= 0.66));                       % true = middle, pred = worst
        WBc = length(find(Yval > 0.66 & fittedcorrected(:,1) <= 0.33));                       % true = worst, pred = best
        WMc = length(find(Yval > 0.66 & (fittedcorrected(:,1) > 0.33 & fittedcorrected(:,1) < 0.66)));      % true = worst, pred = middle
        WWc = length(find(Yval > 0.66 & fittedcorrected(:,1) >= 0.66));                       % true = worst, pred = worst

                
        % store GOF measures 
        GOF.(CV).(F)(1,1) = ModVal.Rsquared.Ordinary(1,1);        % fitted R2
        GOF.(CV).(F)(2,1) = 1-SSE/SSTO;        % fitted R2
        GOF.(CV).(F)(3,1) = RMSEfit;             % RMSE fit
        GOF.(CV).(F)(4,1) = RMSEpred;
        GOF.(CV).(F)(5,1) = RMSEfitcorrected;       % RMSEfitcorr
        GOF.(CV).(F)(6,1) = BB;
        GOF.(CV).(F)(7,1) = BM;
        GOF.(CV).(F)(8,1) = BW;
        GOF.(CV).(F)(9,1) = MB;
        GOF.(CV).(F)(10,1) = MM;
        GOF.(CV).(F)(11,1) = MW;
        GOF.(CV).(F)(12,1) = WB;
        GOF.(CV).(F)(13,1) = WM;
        GOF.(CV).(F)(14,1) = WW;
        GOF.(CV).(F)(15,1) = BBc;
        GOF.(CV).(F)(16,1) = BMc;
        GOF.(CV).(F)(17,1) = BWc;
        GOF.(CV).(F)(18,1) = MBc;
        GOF.(CV).(F)(19,1) = MMc;
        GOF.(CV).(F)(20,1) = MWc;
        GOF.(CV).(F)(21,1) = WBc;
        GOF.(CV).(F)(22,1) = WMc;
        GOF.(CV).(F)(23,1) = WWc;
        GOF.(CV).(F)(24,1) = sum(lastlacval.Corr == 1 & lastlacval.LastLac == 1);
        GOF.(CV).(F)(25,1) = sum(lastlacval.Corr == 1 & lastlacval.LastLac > 1);

        
        T = T+1;
        TT = TT+1;
    end
    clear SSE SSTO RMSEfit RMSEpred PREDBEST PREDWORST WORST BEST iii Xreg Xval Yreg Yval F VN traiindx testindx SSE SSTO x y
end
clear T  resid i ii CV ans ALFIRST RMSE RMSEBEST ModVal SF2 fitted BB BM BW MM MB MW WB WM WW STPWSb STPWSinmodel STPWSstats Yfitted


% Overall performance over the CV
GOF.ALL = array2table((1:25)','VariableNames',{'IDtest'});
GOF.ALL.Name(:,1) = {'R2fit','R2pred','RMSEC','RMSECV','RMSECVcor','BB','BM','BW','MB','MM','MW','WB','WM','WW','BBc','BMc','BWc','MBc','MMc','MWc','WBc','WMc','WWc','FirstCor','HighCor'}';

for ii = 1:10
    CV = sprintf('CV_%d',ii);
    
    AV_R2fit(ii,:) = GOF.(CV)(1,3:15);
    AV_R2pred(ii,:) = GOF.(CV)(2,3:15);
    AV_RMSEfit(ii,:) = GOF.(CV)(3,3:15);
    AV_RMSEpred(ii,:) = GOF.(CV)(4,3:15);
    AV_RMSEpredcorr(ii,:) = GOF.(CV)(5,3:15);
    N_BB(ii,:) = GOF.(CV)(6,3:15);
    N_BM(ii,:) = GOF.(CV)(7,3:15);
    N_BW(ii,:) = GOF.(CV)(8,3:15);
    N_MB(ii,:) = GOF.(CV)(9,3:15);
    N_MM(ii,:) = GOF.(CV)(10,3:15);
    N_MW(ii,:) = GOF.(CV)(11,3:15);
    N_WB(ii,:) = GOF.(CV)(12,3:15);
    N_WM(ii,:) = GOF.(CV)(13,3:15);
    N_WW(ii,:) = GOF.(CV)(14,3:15);
    N_BBc(ii,:) = GOF.(CV)(15,3:15);
    N_BMc(ii,:) = GOF.(CV)(16,3:15);
    N_BWc(ii,:) = GOF.(CV)(17,3:15);
    N_MBc(ii,:) = GOF.(CV)(18,3:15);
    N_MMc(ii,:) = GOF.(CV)(19,3:15);
    N_MWc(ii,:) = GOF.(CV)(20,3:15);
    N_WBc(ii,:) = GOF.(CV)(21,3:15);
    N_WMc(ii,:) = GOF.(CV)(22,3:15);
    N_WWc(ii,:) = GOF.(CV)(23,3:15);
    N_FirstCor(ii,:) = GOF.(CV)(24,3:15);
    N_HighCor(ii,:) = GOF.(CV)(25,3:15);
end

for i = 1:13
    AV_R2pred{find(AV_R2pred{:,i}<0),i} = 0;
end 

GOF.ALL{1,3:15} = nanmean(AV_R2fit{:,:});
GOF.ALL{2,3:15} = nanmean(AV_R2pred{:,:}); 
GOF.ALL{3,3:15} = nanmean(AV_RMSEfit{:,:});
GOF.ALL{4,3:15} = nanmean(AV_RMSEpred{:,:});
GOF.ALL{5,3:15} = nanmean(AV_RMSEpredcorr{:,:});
GOF.ALL{6,3:15} = nansum(N_BB{:,:});
GOF.ALL{7,3:15} = nansum(N_BM{:,:});
GOF.ALL{8,3:15} = nansum(N_BW{:,:});
GOF.ALL{9,3:15} = nansum(N_MB{:,:});
GOF.ALL{10,3:15} = nansum(N_MM{:,:});
GOF.ALL{11,3:15} = nansum(N_MW{:,:});
GOF.ALL{12,3:15} = nansum(N_WB{:,:});
GOF.ALL{13,3:15} = nansum(N_WM{:,:});
GOF.ALL{14,3:15} = nansum(N_WW{:,:});
GOF.ALL{15,3:15} = nansum(N_BBc{:,:});
GOF.ALL{16,3:15} = nansum(N_BMc{:,:});
GOF.ALL{17,3:15} = nansum(N_BWc{:,:});
GOF.ALL{18,3:15} = nansum(N_MBc{:,:});
GOF.ALL{19,3:15} = nansum(N_MMc{:,:});
GOF.ALL{20,3:15} = nansum(N_MWc{:,:});
GOF.ALL{21,3:15} = nansum(N_WBc{:,:});
GOF.ALL{22,3:15} = nansum(N_WMc{:,:});
GOF.ALL{23,3:15} = nansum(N_WWc{:,:});

GOF.ALL{34,3:15} = nansum(N_FirstCor{:,:}); GOF.ALL{34,2} = {'FirstCor'};       
GOF.ALL{35,3:15} = nansum(N_HighCor{:,:});  GOF.ALL{35,2} = {'HighCor'};       

GOF.ALL{24,3:15} = sum(GOF.ALL{6:8,3:15}); GOF.ALL{24,2} = {'B'};
GOF.ALL{25,3:15} = sum(GOF.ALL{9:11,3:15}); GOF.ALL{25,2} = {'M'}; 
GOF.ALL{26,3:15} = sum(GOF.ALL{12:14,3:15}); GOF.ALL{26,2} = {'W'};
GOF.ALL{27,3:15} = sum(GOF.ALL{15:17,3:15}); GOF.ALL{27,2} = {'Bc'};
GOF.ALL{28,3:15} = sum(GOF.ALL{18:20,3:15}); GOF.ALL{28,2} = {'Mc'}; 
GOF.ALL{29,3:15} = sum(GOF.ALL{21:23,3:15}); GOF.ALL{29,2} = {'Wc'};

% true positive and complete opposite %
GOF.ALL{30,3:15} = (GOF.ALL{5,3:15}+GOF.ALL{9,3:15}+GOF.ALL{13,3:15})./sum(GOF.ALL{24:26,3:15}); GOF.ALL{30,2} = {'TRUEPRED'};  GOF.ALL{30,1} = 30;  % 
GOF.ALL{31,3:15} = (GOF.ALL{15,3:15}+GOF.ALL{19,3:15}+GOF.ALL{23,3:15})./sum(GOF.ALL{27:29,3:15}); GOF.ALL{31,2} = {'TRUEPREDc'};  GOF.ALL{31,1} = 31;  % 
GOF.ALL{32,3:15} = (GOF.ALL{8,3:15}+GOF.ALL{12,3:15})./sum(GOF.ALL{24:26,3:15}); GOF.ALL{32,2} = {'OPPOSITE'};  GOF.ALL{32,1} = 32;
GOF.ALL{33,3:15} = (GOF.ALL{17,3:15}+GOF.ALL{21,3:15})./sum(GOF.ALL{24:26,3:15}); GOF.ALL{33,2} = {'OPPOSITEc'};  GOF.ALL{33,1} = 33;

% correctly predicted and in first lactation
GOF.ALL{36,3:15} = GOF.ALL{34,3:15}./(GOF.ALL{15,3:15}+GOF.ALL{19,3:15}+GOF.ALL{23,3:15}); GOF.ALL{36,2} = {'Correct1st'};  % percentage correctly predicted in first lactation
GOF.ALL{37,3:15} = GOF.ALL{35,3:15}./(GOF.ALL{15,3:15}+GOF.ALL{19,3:15}+GOF.ALL{23,3:15}); GOF.ALL{37,2} = {'CorrectHigher'};  % percentage correctly predicted in higher lactation

% total number of 1st lactation = last lactation cows in test sets
for ii = 1:10
    T = 1;
    for i = Farms
        testindx = find(ismember(find(X.FarmID == i),test(:,ii))==1);   % select indices of testset

        y = sortrows(Y(Y.FarmID == i,:),2);     % select y
        
        Nfirstlact = y(testindx,:);         % ytraining
        Testsetlact1(T,ii) = length(find(Nfirstlact.LastLac == 1));
        Testsetlact2(T,ii) = length(find(Nfirstlact.LastLac > 1));
        T = T+1;
    end
end
percfirsttestset = Testsetlact1./(Testsetlact1+Testsetlact2);
percfirsttestset(:,11) = mean(percfirsttestset(:,1:10),2)*100;

Testsetlact1(:,11) = sum(Testsetlact1,2);
correct1s = GOF.ALL{34,3:end}'./Testsetlact1(:,11);

clear N_TH N_TH_PH N_TL N_TL_PH i  ii fitted CV
clear N_BB N_BM N_BW N_MM N_MW N_MB N_WW N_WB N_WM N_BBc N_BMc N_BWc N_MMc N_MWc N_MBc N_WWc N_WBc N_WMc
clear WBc WMc WWc MBc MMc MWc BBc BMc BWc B A TrainFit  RMSEfitcorr residcorrected fittedcorrected AV_R2fit AV_R2pred AV_RMSEpred AV_RMSEpredcorr AV_RMSEfit RMSEfitcorrected



