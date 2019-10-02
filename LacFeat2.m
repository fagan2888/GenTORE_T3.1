function [OUT,OUT2,OUT3,MINOR,MAJOR] = LacFeat2(DATA,D)
% this function calculates SF
%       WOOD after iterations
%           A, B, C
%           RMSE
%           Peak yield
%           DIM peak
%           Peak/DIM peak
%           Perturbations characteristics
%       WOOD first model
%           persistency
%       305 D yield % ct average herd
%
% DATA contains at least CowId Lac TDMY 
%
% STEP 1: calculate features for each cow/lac in DATA
% STEP 2: calculate herd averages
% STEP 3: correct features for herd averages

% DATA = DAY_LT.VanDenHaute; D = 0;

OUT = unique([DATA.CowId DATA.Lac],'rows');
OUT = array2table(OUT,'VariableNames',{'CowId','Lac'});
OUT.Ndays(:,1) = NaN;           % number of the days the lactation lasts
OUT.Y305D(:,1) = NaN;           % total yield of at most 305 days 
OUT.Y305Dc(:,1) = NaN;          % yield further calculated to 305 days
OUT.RMSEi(:,1) = NaN;           % initial RMSE
OUT.PeakYi(:,1) = NaN;          % peak based on Wood initial
OUT.PeakDi(:,1) = NaN;          % peak based on Wood initial (first within 95%)
OUT.Ratei(:,1) = NaN;           % rate of increase based on Wood initial
OUT.Persi(:,1) = NaN;           % persistency based on Wood initial
OUT.RMSEc(:,1) = NaN;           % RMSE after iterations
OUT.PeakYc(:,1) = NaN;          % peak based on Wood after iterations
OUT.PeakDc(:,1) = NaN;          % peak DIM based on Wood (first within 95%) after iterations
OUT.Ratec(:,1) = NaN;           % rate of increase based on Wood after iterations
OUT.Persc(:,1) = NaN;           % persistency based on Wood initial
OUT.pLOSSc(:,1) = NaN;          % sum of the strong negative yield residuals
OUT.pNEGc(:,1) = NaN;           % ratio positive vs negative residuals
OUT.NPertMJ(:,1) = NaN;         % number of perturbations
OUT.RecDaysMJ(:,1) = NaN;       % average recovery time major perturbations
OUT.RecLossMJ(:,1) = NaN;       % average 'loss' associated with perturbation recovery
OUT.DevLossMJ(:,1) = NaN;       % average loss associated with perturbation development
OUT.DevDaysMJ(:,1) = NaN;       % average development time major perturbations
OUT.AvMinMJ(:,1) = NaN;         % average minimum of detected losses
OUT.NPertMI(:,1) = NaN;         % number of perturbations
OUT.RecDaysMI(:,1) = NaN;       % average N days needed for recovery
OUT.RecLossMI(:,1) = NaN;       % average 'loss' associated with perturbation recovery
OUT.DevLossMI(:,1) = NaN;       % average loss associated with perturbation development
OUT.DevDaysMI(:,1) = NaN;       % average development time minor perturbations
OUT.AvMinMI(:,1) = NaN;         % average minimum of detected losses

OUT.N5dayNeg(:,1) = NaN;        % N phases of > 5 days negative residuals and at least once below 85 of model
OUT.N10dayNeg(:,1) = NaN;       % N phases of > 10 days negative residuals and at least once below 85% of model
OUT.N1dayNes(:,1) = NaN;        % percentage of days < 0.85% Wood in first 305days
OUT.ARes3per(:,1) = NaN;        % Average "loss" of (3) largest perturbations when production compared to WOOD model  
OUT.maxRes(:,1) = NaN;         % Absolute max negative residual within period of at least 5 negatives
OUT.Autocorr(:,1) = NaN;       % 'Autocorrelation' expressed as number of sign changes in the residuals
OUT.MeanRes(:,1) = NaN;        % Mean of the residuals Wood
OUT.AbsMeanRes(:,1) = NaN;     % Mean of the absolute value of the residuals Wood
OUT.VarRes(:,1) = NaN;         % Variance of the residuals Wood
OUT.RatioRMSE(:,1) = NaN;       % ratio between two RMSE

% Distinct for timing-wise SF
OUT.Y100D(:,1) = NaN;           % Total yield of at most 305 days 
OUT.RMSEi100(:,1) = NaN;        % fill in RSME initial 
OUT.PeakYi100(:,1) = NaN;       % peak based on Wood
OUT.PeakDi100(:,1) = NaN;       % peak based on Wood (first within 95%)
OUT.Ratei100(:,1) = NaN;        % rate of increase
OUT.RMSEc100(:,1) = NaN;        % RMSE of Wood model after correction in first 100 days
OUT.RatioRMSE100(:,1) = NaN;    % Ratio of RMSEs
OUT.pLOSSc100(:,1) = NaN;       % sum of the strong negative yield residuals
OUT.pNEGc100(:,1) = NaN;        % ratio negative and positive residuals
OUT.maxRes100(:,1) = NaN;       % minimal residual
OUT.MeanRes100(:,1) = NaN;      % Mean of the residuals Wood
OUT.AbsMeanRes100(:,1) = NaN;   % Mean of the absolute value of the residuals Wood
OUT.VarRes100(:,1) = NaN;       % Variance of the residuals Wood   
OUT.Autocorr100(:,1) = NaN;     % 'Autocorrelation' expressed as number of sign changes in the residuals vs no of res
OUT.NPertMJ100(:,1) = NaN;      % number of perturbations
OUT.RecDaysMJ100(:,1) = NaN;    % average recovery time major perturbations
OUT.RecLossMJ100(:,1) = NaN;    % average 'loss' associated with perturbation recovery
OUT.DevLossMJ100(:,1) = NaN;    % average loss associated with perturbation development
OUT.DevDaysMJ100(:,1) = NaN;    % average development time major perturbations
OUT.AvMinMJ100(:,1) = NaN;      % average minimum of detected losses
OUT.NPertMI100(:,1) = NaN;      % number of perturbations
OUT.RecDaysMI100(:,1) = NaN;    % average N days needed for recovery
OUT.RecLossMI100(:,1) = NaN;    % average 'loss' associated with perturbation recovery
OUT.DevLossMI100(:,1) = NaN;    % average loss associated with perturbation development
OUT.DevDaysMI100(:,1) = NaN;    % average development time minor perturbations
OUT.AvMinMI100(:,1) = NaN;      % average minimum of detected losses
OUT.NPertALL100(:,1) = NaN;     % number of perturbations
OUT.RecDaysALL100(:,1) = NaN;   % average recovery time all perturbations
OUT.RecLossALL100(:,1) = NaN;   % average 'loss' associated with perturbation recovery
OUT.DevLossALL100(:,1) = NaN;   % average loss associated with perturbation development
OUT.DevDaysALL100(:,1) = NaN;   % average development time all perturbations
OUT.AvMinALL100(:,1) = NaN;     % average min of all perturbations

% prepare perturbation output
MINOR = array2table(zeros(10000,17),'VariableNames',{'CowId','Lac','Npert','TMY','DIM','IDXstart','DIMstart','IDXend','DIMend','Length','MinLos','IDXmin','DIMmin','RECO','Loss','DEVE','DLOSS'});
MAJOR = array2table(zeros(10000,17),'VariableNames',{'CowId','Lac','Npert','TMY','DIM','IDXstart','DIMstart','IDXend','DIMend','Length','MinLos','IDXmin','DIMmin','RECO','Loss','DEVE','DLOSS'});
M1 = 1; M2 = 1;

options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',5000,'Display','off');    % options
% close all
Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t); % Wood model

for i = 1:length(OUT.Lac)
    
    ind = find(DATA.CowId == OUT.CowId(i,1) & DATA.Lac == OUT.Lac(i,1) & datenum(DATA.Date)-datenum(DATA.Calving) <= 305 & isnan(DATA.TDMY)==0);
    idx = find(datenum(DATA.Date(ind)) - datenum(DATA.Calving(ind))>=0);
    ind = ind(idx); clear idx
    
    if length(ind) > 200
        TMY = DATA.TDMY(ind);                       % TDMY
        DIM = datenum(DATA.Date(ind))- datenum(DATA.Calving(ind)); % DIM

            OUT.Ndays(i,1) = DIM(end);              % Number of the days the lactation lasts
            OUT.Y305D(i,1) = sum(TMY);              % Total yield of at most 305 days 
            OUT.Y305Dc(i,1) = (sum(TMY)/length(DIM))*305; % 305 day yield corrected

        % if DELAVAL, first smooth data (moving average; span = 3)
        if D == 1
            TMY = smooth(TMY,5);                    % smoothed
        end 

        % Wood before any iterations

        pWOOD = lsqcurvefit(Wood,[15 0.20 0.004],DIM,TMY,[0 0 0],[20 1 0.01],options);
        fWOOD = Wood(pWOOD,DIM);                    % function values of Wood
        RMSE = sqrt(mean((TMY-fWOOD).^2));          % RMSE of Wood model
        rWOOD = TMY - fWOOD;                        % residuals of Wood model
            OUT.RMSEi(i,1) = RMSE;                  % fill in RSME initial 
            OUT.PeakYi(i,1) = max(Wood(pWOOD,DIM));      % peak based on Wood
            OUT.PeakDi(i,1) = DIM(find(Wood(pWOOD,DIM) > 0.99* max(Wood(pWOOD,DIM)),1,'first'));      % peak based on Wood (first within 95%)
            OUT.Ratei(i,1) = OUT.PeakYi(i,1)./OUT.PeakDi(i,1);  % rate of increase
            OUT.Persi(i,1) = (Wood(pWOOD,200)-Wood(pWOOD,180))/20; 


% % % %                         figure(i); hold on;
% % % %                         plot(DIM, TMY, 'o-','LineWidth',2);
% % % %                         plot(DIM, fWOOD,'-','LineWidth',2);
% % % %                         plot(OUT.PeakDi(i,1),OUT.PeakYi(i,1),'*','MarkerSize',10,'LineWidth',2);
% % % %                         plot([180 200],[Wood(pWOOD,180) Wood(pWOOD,200)],'--','LineWidth',2);


        % Wood after iterations
        DIM2 = DIM;                                     % prepare DIM
        TMY2 = TMY;                                     % prepare TMY
        fWOOD2 = fWOOD;                                 % prepare fWOOD
        p = pWOOD;                                      % prepare parameters
        RMSE2 = 0;                                      % prepare RMSE     
        n=0;                                            % Counter for the number of iterations
        while abs(RMSE-RMSE2) > 0.1  && n <= 20         % while the difference between the model and the previous model is larger than 0.1 (??)
            n=n+1;                                      % and also, at most 20 iterations
            DIM2 = DIM2(TMY2>0.95*fWOOD2);              % DIM of all milkings above 0.85 of the prediction
            TMY2 = TMY2(TMY2>0.95*fWOOD2);              % TMY of all milkings above 0.85 of the prediction

            if length(DIM2)>4
                p = lsqcurvefit(Wood,p,DIM2,TMY2,[0 0 0],[20 1 0.01],options);  % refit model
                fWOOD2 = Wood(p,DIM2);                  % calculate model values

                RMSE2 = RMSE;                           % Set the RMSE to the previous ones
                RMSE = sqrt(mean((TMY2-fWOOD2).^2));    % RMSE
            else
                n = 100;
            end
        end


        % lactation curve characteristics based on Wood after iterations (c)
            OUT.RMSEc(i,1) = sqrt(mean((TMY-Wood(p,DIM)).^2));      % fill in RSME after iterations (corr)
            OUT.PeakYc(i,1) = max(Wood(p,DIM));         % peak based on Wood
            OUT.PeakDc(i,1) = DIM(find(Wood(p,DIM) > 0.99* max(Wood(p,DIM)),1,'first'));      % peak based on Wood (first within 95%)
            OUT.Ratec(i,1) = OUT.PeakYc(i,1)./OUT.PeakDc(i,1);      % rate of increase
            OUT.Persc(i,1) = (Wood(p,200)-Wood(p,180))/20;          % persistency based on Wood initial


% % % %                         plot(DIM, Wood(p,DIM),'.-','LineWidth',2);
% % % %                         plot(OUT.PeakDc(i,1),OUT.PeakYc(i,1),'*','MarkerSize',10,'LineWidth',2);
% % % %                         plot([180 200],[Wood(p,180) Wood(p,200)],'--','LineWidth',2);


        % perturbation characteristics - based on Wood after iterations (c)
        idx = find(TMY<0.85*Wood(p,DIM));
            OUT.pLOSSc(i,1) = sum(TMY(idx)-Wood(p,DIM(idx)));       % sum of the strong negative yield residuals

% % % %                         plot(DIM(idx), TMY(idx),'ro','LineWidth',2); % the residuals below 90% of the prediction


        idx1 = find(TMY<Wood(p,DIM));                               % all negative residuals
        idx2 = find(TMY>Wood(p,DIM));                               % all positive residuals

            OUT.pNEGc(i,1) = abs(sum(TMY(idx2)-Wood(p,DIM(idx2)))./sum(TMY(idx1)-Wood(p,DIM(idx1))))*100; % ratio negative and positive residuals
            
        
        % recovery rate - based on Wood after iterations - MAJOR
        pert1 = IsoPert(DIM,TMY,Wood(p,DIM),20,10,40,2,2);          % call perturbation function - A=2 (procentual loss) and P=2, bended regression
                        
            OUT.NPertMJ(i,1) = max(pert1.Npert);                    % number of perturbations
            if OUT.NPertMJ(i,1) ~= 0
                OUT.RecDaysMJ(i,1) = nanmean(pert1.RECO);           % average recovery time major perturbations
                OUT.RecLossMJ(i,1) = nanmean(pert1.Loss);           % average 'loss' associated with perturbation recovery
                OUT.DevLossMJ(i,1) = nanmean(pert1.DLOSS);          % average loss associated with perturbation development
                OUT.DevDaysMJ(i,1) = nanmean(pert1.DEVE);           % average development time major perturbations
                OUT.AvMinMJ(i,1) = nanmean(pert1.MinLos);           % average minimum of detected losses
            end
               
        % recovery rate - based on Wood after iterations - MINOR 
        pert2 = IsoPert(DIM,TMY,Wood(p,DIM),10,5,30,2,2);           % call perturbation function
        idx = find(ismember(pert2.IDXstart,pert1.IDXstart)==1 & ismember(pert2.IDXend,pert1.IDXend)==1);
        pert2(idx,:) = [];                                          % delete major perturbations
        if isempty(pert2)== 0; pert2.Npert(:,1) = 1:length(pert2.Npert); end % adjust N perturbations
                        
            if isempty(pert2) == 0
                OUT.NPertMI(i,1) = max(pert2.Npert);                % number of perturbations
                OUT.RecDaysMI(i,1) = nanmean(pert2.RECO);           % average N days needed for recovery
                OUT.RecLossMI(i,1) = nanmean(pert2.Loss);           % average 'loss' associated with perturbation recovery
                OUT.DevLossMI(i,1) = nanmean(pert2.DLOSS);          % average loss associated with perturbation development
                OUT.DevDaysMI(i,1) = nanmean(pert2.DEVE);           % average development time minor perturbations
                OUT.AvMinMI(i,1) = nanmean(pert2.MinLos);           % average minimum of detected losses
            else 
                OUT.NPertMI(i,1) = 0;                % number of perturbations
            end
            
        % further characterization of the negative residuals    
        fWOOD2 = Wood(p,DIM);                                       % Wood model after asymmetric least squares
        rWOOD2 = TMY-fWOOD;                                         % Residuals of the Wood model
        idx = find(TMY<0.85*fWOOD2);                                % residuals smaller than 0.85 of wood
        NegRes = (rWOOD2<0);                                        % 1 for all negative Wood residuals
        [~, idx2] = Pattern(NegRes,[1 1 1 1 1]);                    % indices for which 5 consecutive measurements the residual is negative
        Negs = zeros(length(idx2),6);
        for j = 1:length(idx2)                                      % for all the periods in which 5 conse negative residuals
            idx3 = find(DIM>DIM(idx2(j)) & rWOOD2 > 0,1,'first');   % endpoint DIM
            if isempty(idx3); idx3 = length(rWOOD2);end             % no recovery
            if idx3 ~= length(rWOOD2)                               % if recovery within the lactation (and 305 days)
                Negs(j,1:2) = [DIM(idx2(j)) DIM(idx3)-1 ];          % DIM start and DIM end
            else
                Negs(j,1:2) = [DIM(idx2(j)) DIM(idx3)];             % DIM start and DIM end
            end
            Negs(j,3) = max(ismember(idx2(j):idx3,idx));            
            Negs(j,4) = DIM(idx3)-1-DIM(idx2(j));
            Negs(j,5) = sum(TMY(idx2(j):idx3)-fWOOD2(idx2(j):idx3)); % 
        end
        
            OUT.N5dayNeg(i,1) = nansum(Negs(:,3));                  % N phases of > 5 days negative residuals and at least once below 85% model value
            OUT.N10dayNeg(i,1) = sum(Negs(Negs(:,4)>=10,3));        % N phases of > 10 days negative residuals and at least once below 85% model value
            OUT.N1dayNes(i,1) = (length(idx)/length(DIM))*100;      % percentage of residuals below 85% of wood compared to N of measurments
        
        Negs = sortrows(Negs, 5,'ascend');                          % sort by 'size' of perturbation in terms of total loss compared WOOD
            OUT.ARes3per(i,1) = sum(Negs(1:min(3,j),5))/sum(Negs(1:min(3,j),4));    % Average "loss" of (3) largest perturbations when production compared to WOOD model  

            OUT.maxRes(i,1) = min(rWOOD2(4:end));                    % minimal residual

            OUT.MeanRes(i,1) = nanmean(rWOOD2);                      % Mean of the residuals Wood
            OUT.AbsMeanRes(i,1) = nanmean(abs(rWOOD2));              % Mean of the absolute value of the residuals Wood
            OUT.VarRes(i,1) = nanvar(rWOOD2);                        % Variance of the residuals Wood   
            
            OUT.Autocorr(i,1) = (1-length(find(diff(sign(rWOOD2))>1))/length(rWOOD2))*100; % 'Autocorrelation' expressed as number of sign changes in the residuals vs no of res
            
        % difference between model iterated and model corrected
            OUT.RatioRMSE(i,1) = (OUT.RMSEc(i,1)/OUT.RMSEi(i,1))*100; % ratio between two RMSE
       
       
        % Store and output perturbations per lactation
        L1 = length(pert1.Npert);                                    % number of major perturbations
        L2 = length(pert2.Npert);                                    % number of minor perturbations
        MAJOR(M1:M1+L1-1,3:17) = pert1;                              % fill in major data
        MINOR(M2:M2+L2-1,3:17) = pert2;                              % fill in minor data
        MAJOR.CowId(M1:M1+L1-1,1) = OUT.CowId(i);                    % fill in cowId
        MAJOR.Lac(M1:M1+L1-1,1) = OUT.Lac(i);                        % fill in Lac
        MINOR.CowId(M2:M2+L2-1,1) = OUT.CowId(i);                    % fill in cowId
        MINOR.Lac(M2:M2+L2-1,1) = OUT.Lac(i);                        % fill in Lac
        M1 = M1+L1; M2 = M2+L2;                                      % increase tellers
            
       
        % Timing-based sensor features, to calculate on the first 100 days in lactation     
            % fit a new model which only uses DIM and TMY of first 100 days
        ind = find(DATA.CowId == OUT.CowId(i,1) & DATA.Lac == OUT.Lac(i,1) & datenum(DATA.Date)-datenum(DATA.Calving) <= 100 & isnan(DATA.TDMY)==0);
        TMY = DATA.TDMY(ind);                                       % TDMY
        DIM = datenum(DATA.Date(ind))- datenum(DATA.Calving(ind));  % DIM

            OUT.Y100D(i,1) = sum(TMY);                              % Total yield of at most 305 days 
        
        % if DELAVAL, first smooth data (moving average; span = 3)
        if D == 1
            TMY = smooth(TMY,5);                                    % smoothed
        end
        
        pWOOD1 = lsqcurvefit(Wood,[15 0.20 0.004],DIM,TMY,[0 0 0],[25 2 0.05],options);
        fWOOD1 = Wood(pWOOD1,DIM);                                  % function values of Wood
        RMSE1 = sqrt(mean((TMY-fWOOD1).^2));                        % RMSE of Wood model
        rWOOD1 = TMY - fWOOD1;                                      % residuals of Wood model
            OUT.RMSEi100(i,1) = RMSE1;                              % fill in RSME initial 
            OUT.PeakYi100(i,1) = max(fWOOD1);                       % peak based on Wood
            OUT.PeakDi100(i,1) = DIM(find(fWOOD1 > 0.99* max(fWOOD1),1,'first')); % peak based on Wood (first within 95%)
            OUT.Ratei100(i,1) = OUT.PeakYi100(i,1)./OUT.PeakDi100(i,1); % rate of increase

        % Wood after iterations
        DIM3 = DIM;                                         % prepare DIM
        TMY3 = TMY;                                         % prepare TMY
        fWOOD3 = fWOOD1;                                    % prepare fWOOD
        p = pWOOD1;                                         % prepare parameters
        RMSE3 = 0;                                          % prepare RMSE     
        n=0;                                                % Counter for the number of iterations
        while abs(RMSE1-RMSE3) > 0.1  && n <= 20            % while the difference between the model and the previous model is larger than 0.1 (??)
            n=n+1;                                          % and also, at most 20 iterations
            DIM3 = DIM3(TMY3>0.95*fWOOD3);                  % DIM of all milkings above 0.85 of the prediction
            TMY3 = TMY3(TMY3>0.95*fWOOD3);                  % TMY of all milkings above 0.85 of the prediction

            if length(DIM3)>4
                p = lsqcurvefit(Wood,p,DIM3,TMY3,[0 0 0],[20 1 0.01],options);  % refit model
                fWOOD3 = Wood(p,DIM3);                      % calculate model values

                RMSE3 = RMSE1;                              % Set the RMSE to the previous ones
                RMSE1 = sqrt(mean((TMY3-fWOOD3).^2));       % RMSE
            else
                n = 100;
            end
        end
        
        fWOOD3 = Wood(p,DIM);                               % Wood of iterated
        rWOOD3 = TMY-fWOOD3;                                % residuals of Wood iterated
            
            OUT.RMSEc100(i,1) = sqrt(mean((TMY-fWOOD3).^2));% RMSE of Wood model after correction in first 100 days
            OUT.RatioRMSE100(i,1) = (OUT.RMSEc100(i,1)/OUT.RMSEi100(i,1))*100;   % Ratio of RMSEs

        idx = find(TMY<0.85*fWOOD3);
            OUT.pLOSSc100(i,1) = sum(TMY(idx)-fWOOD3(idx)); % sum of the strong negative yield residuals

        idx1 = find(TMY<fWOOD3);                            % all negative residuals
        idx2 = find(TMY>fWOOD3);                            % all positive residuals
            OUT.pNEGc100(i,1) = abs(sum(TMY(idx2)-Wood(p,DIM(idx2)))./sum(TMY(idx1)-Wood(p,DIM(idx1))))*100; % ratio negative and positive residuals
   
            OUT.maxRes100(i,1) = min(rWOOD3(4:end));        % minimal residual

            OUT.MeanRes100(i,1) = nanmean(rWOOD3);          % Mean of the residuals Wood
            OUT.AbsMeanRes100(i,1) = nanmean(abs(rWOOD3));  % Mean of the absolute value of the residuals Wood
            OUT.VarRes100(i,1) = nanvar(rWOOD3);            % Variance of the residuals Wood   
            
            OUT.Autocorr100(i,1) = (1-length(find(diff(sign(rWOOD3))>1))/length(rWOOD3))*100; % 'Autocorrelation' expressed as number of sign changes in the residuals vs no of res

        % select and characterize the minor and major perturbations in first 100 days
        % recovery rate - based on Wood after iterations - MAJOR
        pert3 = IsoPert(DIM,TMY,fWOOD3,20,5,40,2,2);       % call perturbation function - A=2 (procentual loss) and P=2, bended regression
                        
            OUT.NPertMJ100(i,1) = max(pert3.Npert);                    % number of perturbations
            if OUT.NPertMJ100(i,1) ~= 0
                OUT.RecDaysMJ100(i,1) = nanmean(pert3.RECO);           % average recovery time major perturbations
                OUT.RecLossMJ100(i,1) = nanmean(pert3.Loss);           % average 'loss' associated with perturbation recovery
                OUT.DevLossMJ100(i,1) = nanmean(pert3.DLOSS);          % average loss associated with perturbation development
                OUT.DevDaysMJ100(i,1) = nanmean(pert3.DEVE);           % average development time major perturbations
                OUT.AvMinMJ100(i,1) = nanmean(pert3.MinLos);           % average minimum of detected losses
            end
            
        % recovery rate - based on Wood after iterations - MINOR 
        pert4 = IsoPert(DIM,TMY,Wood(p,DIM),10,5,30,2,2);              % call perturbation function
        idx = find(ismember(pert4.IDXstart,pert3.IDXstart)==1 & ismember(pert4.IDXend,pert3.IDXend)==1);
        pert4(idx,:) = [];                                             % delete major perturbations
        if isempty(pert4)== 0; pert4.Npert(:,1) = 1:length(pert4.Npert); end % adjust N perturbations
                        
            if isempty(pert4) == 0
                OUT.NPertMI100(i,1) = max(pert4.Npert);                % number of perturbations
                OUT.RecDaysMI100(i,1) = nanmean(pert4.RECO);           % average N days needed for recovery
                OUT.RecLossMI100(i,1) = nanmean(pert4.Loss);           % average 'loss' associated with perturbation recovery
                OUT.DevLossMI100(i,1) = nanmean(pert4.DLOSS);          % average loss associated with perturbation development
                OUT.DevDaysMI100(i,1) = nanmean(pert4.DEVE);           % average development time minor perturbations
                OUT.AvMinMI100(i,1) = nanmean(pert4.MinLos);           % average minimum of detected losses
            else 
                OUT.NPertMI100(i,1) = 0;                               % number of perturbations
            end
       % all pertrbations combined
       pert5 = IsoPert(DIM,TMY,Wood(p,DIM),10,5,30,2,2);               % call perturbation function

            OUT.NPertALL100(i,1) = max(pert5.Npert);                   % number of perturbations
            if OUT.NPertALL100(i,1) ~= 0
                OUT.RecDaysALL100(i,1) = nanmean(pert5.RECO);              % average recovery time major perturbations
                OUT.RecLossALL100(i,1) = nanmean(pert5.Loss);              % average 'loss' associated with perturbation recovery
                OUT.DevLossALL100(i,1) = nanmean(pert5.DLOSS);             % average loss associated with perturbation development
                OUT.DevDaysALL100(i,1) = nanmean(pert5.DEVE);              % average development time major perturbations
                OUT.AvMinALL100(i,1) = nanmean(pert5.MinLos);              % average minimum of detected losses
            end
            
    elseif length(ind)>100              % only calculate chars for first 100 days
        
        ind = find(DATA.CowId == OUT.CowId(i,1) & DATA.Lac == OUT.Lac(i,1) & datenum(DATA.Date)-datenum(DATA.Calving) <= 100 & isnan(DATA.TDMY)==0);
        TMY = DATA.TDMY(ind);                                       % TDMY
        DIM = datenum(DATA.Date(ind))- datenum(DATA.Calving(ind));  % DIM

            OUT.Y100D(i,1) = sum(TMY);                              % Total yield of at most 305 days 
        
        % if DELAVAL, first smooth data (moving average; span = 3)
        if D == 1
            TMY = smooth(TMY,5);                                    % smoothed
        end
        
        pWOOD1 = lsqcurvefit(Wood,[15 0.20 0.004],DIM,TMY,[0 0 0],[25 2 0.05],options);
        fWOOD1 = Wood(pWOOD1,DIM);                                  % function values of Wood
        RMSE1 = sqrt(mean((TMY-fWOOD1).^2));                        % RMSE of Wood model
        rWOOD1 = TMY - fWOOD1;                                      % residuals of Wood model
            OUT.RMSEi100(i,1) = RMSE1;                              % fill in RSME initial 
            OUT.PeakYi100(i,1) = max(fWOOD1);                       % peak based on Wood
            OUT.PeakDi100(i,1) = DIM(find(fWOOD1 > 0.99* max(fWOOD1),1,'first')); % peak based on Wood (first within 95%)
            OUT.Ratei100(i,1) = OUT.PeakYi100(i,1)./OUT.PeakDi100(i,1); % rate of increase

        % Wood after iterations
        DIM3 = DIM;                                         % prepare DIM
        TMY3 = TMY;                                         % prepare TMY
        fWOOD3 = fWOOD1;                                    % prepare fWOOD
        p = pWOOD1;                                         % prepare parameters
        RMSE3 = 0;                                          % prepare RMSE     
        n=0;                                                % Counter for the number of iterations
        while abs(RMSE1-RMSE3) > 0.1  && n <= 20            % while the difference between the model and the previous model is larger than 0.1 (??)
            n=n+1;                                          % and also, at most 20 iterations
            DIM3 = DIM3(TMY3>0.95*fWOOD3);                  % DIM of all milkings above 0.85 of the prediction
            TMY3 = TMY3(TMY3>0.95*fWOOD3);                  % TMY of all milkings above 0.85 of the prediction

            if length(DIM3)>4
                p = lsqcurvefit(Wood,p,DIM3,TMY3,[0 0 0],[20 1 0.01],options);  % refit model
                fWOOD3 = Wood(p,DIM3);                      % calculate model values

                RMSE3 = RMSE1;                              % Set the RMSE to the previous ones
                RMSE1 = sqrt(mean((TMY3-fWOOD3).^2));       % RMSE
            else
                n = 100;
            end
        end
        
        fWOOD3 = Wood(p,DIM);                               % Wood of iterated
        rWOOD3 = TMY-fWOOD3;                                % residuals of Wood iterated
            
            OUT.RMSEc100(i,1) = sqrt(mean((TMY-fWOOD3).^2));% RMSE of Wood model after correction in first 100 days
        
        idx = find(TMY<0.85*fWOOD3);
            OUT.pLOSSc100(i,1) = sum(TMY(idx)-fWOOD3(idx)); % sum of the strong negative yield residuals

        idx1 = find(TMY<fWOOD3);                            % all negative residuals
        idx2 = find(TMY>fWOOD3);                            % all positive residuals
            OUT.pNEGc100(i,1) = abs(sum(TMY(idx2)-Wood(p,DIM(idx2)))./sum(TMY(idx1)-Wood(p,DIM(idx1))))*100; % ratio negative and positive residuals
   
            OUT.maxRes100(i,1) = min(rWOOD3(4:end));        % minimal residual

            OUT.MeanRes100(i,1) = nanmean(rWOOD3);          % Mean of the residuals Wood
            OUT.AbsMeanRes100(i,1) = nanmean(abs(rWOOD3));  % Mean of the absolute value of the residuals Wood
            OUT.VarRes100(i,1) = nanvar(rWOOD3);            % Variance of the residuals Wood   
            
            OUT.Autocorr100(i,1) = (1-length(find(diff(sign(rWOOD3))>1))/length(rWOOD3))*100; % 'Autocorrelation' expressed as number of sign changes in the residuals vs no of res

        % select and characterize the minor and major perturbations in first 100 days
        % recovery rate - based on Wood after iterations - MAJOR
        pert3 = IsoPert(DIM,TMY,fWOOD3,20,5,40,2,2);       % call perturbation function - A=2 (procentual loss) and P=2, bended regression
                        
            OUT.NPertMJ100(i,1) = max(pert3.Npert);                    % number of perturbations
            if OUT.NPertMJ100(i,1) ~=0
                OUT.RecDaysMJ100(i,1) = nanmean(pert3.RECO);           % average recovery time major perturbations
                OUT.RecLossMJ100(i,1) = nanmean(pert3.Loss);           % average 'loss' associated with perturbation recovery
                OUT.DevLossMJ100(i,1) = nanmean(pert3.DLOSS);          % average loss associated with perturbation development
                OUT.DevDaysMJ100(i,1) = nanmean(pert3.DEVE);           % average development time major perturbations
                OUT.AvMinMJ100(i,1) = nanmean(pert3.MinLos);           % average minimum of detected losses
            end
            
        % recovery rate - based on Wood after iterations - MINOR 
        pert4 = IsoPert(DIM,TMY,Wood(p,DIM),10,5,30,2,2);              % call perturbation function
        idx = find(ismember(pert4.IDXstart,pert3.IDXstart)==1 & ismember(pert4.IDXend,pert3.IDXend)==1);
        pert4(idx,:) = [];                                             % delete major perturbations
        if isempty(pert4)== 0; pert4.Npert(:,1) = 1:length(pert4.Npert); end % adjust N perturbations
                        
            
            if isempty(pert4)== 0
                OUT.NPertMI100(i,1) = max(pert4.Npert);                % number of perturbations
                OUT.RecDaysMI100(i,1) = nanmean(pert4.RECO);           % average N days needed for recovery
                OUT.RecLossMI100(i,1) = nanmean(pert4.Loss);           % average 'loss' associated with perturbation recovery
                OUT.DevLossMI100(i,1) = nanmean(pert4.DLOSS);          % average loss associated with perturbation development
                OUT.DevDaysMI100(i,1) = nanmean(pert4.DEVE);           % average development time minor perturbations
                OUT.AvMinMI100(i,1) = nanmean(pert4.MinLos);           % average minimum of detected losses
            else
                OUT.NPertMI100(i,1) = 0;                                % number of perturbations
            end
            
       % all pertrbations combined
       pert5 = IsoPert(DIM,TMY,Wood(p,DIM),10,5,30,2,2);               % call perturbation function

            OUT.NPertALL100(i,1) = max(pert5.Npert);                   % number of perturbations
            if OUT.NPertALL100(i,1) ~=0
                OUT.RecDaysALL100(i,1) = nanmean(pert5.RECO);          % average recovery time major perturbations
                OUT.RecLossALL100(i,1) = nanmean(pert5.Loss);          % average 'loss' associated with perturbation recovery
                OUT.DevLossALL100(i,1) = nanmean(pert5.DLOSS);         % average loss associated with perturbation development
                OUT.DevDaysALL100(i,1) = nanmean(pert5.DEVE);          % average development time major perturbations
                OUT.AvMinALL100(i,1) = nanmean(pert5.MinLos);          % average minimum of detected losses       
            end
    end
end
clear ind idx n idx1 idx2 DIM DIM2 ans D i p pWOOD RMSE RMSE2 TMY TMY2 fWOOD fWOOD2 rWOOD

% delete surplus rows
MINOR(MINOR.CowId == 0,:) = [];
MAJOR(MAJOR.CowId == 0,:) = [];
%% STEP 2: herd characteristics

OUT2 = unique(OUT.Lac); OUT2(:,2:length(OUT.Properties.VariableNames)) = 0;               % OUT2 will contain per lactation chars
STD = unique(OUT.Lac); STD(:,2:length(OUT.Properties.VariableNames)) = 0;                 % STD will contain per lactation chars std
for i = 1:max(OUT2(:,1))            
    OUT2(i,2) = sum(OUT.Lac == i);                      % number of lactations in ith lac for this herd
    OUT2(i,3:length(OUT.Properties.VariableNames)) = nanmean(OUT{OUT.Lac == i,3:length(OUT.Properties.VariableNames)},1);   % averages per lactation number
    STD(i,3:length(OUT.Properties.VariableNames)) = nanstd(OUT{OUT.Lac == i,3:length(OUT.Properties.VariableNames)},1);     % std per lactation number
end
OUT2 = array2table(OUT2,'VariableNames',OUT.Properties.VariableNames); % array to table
OUT2.Properties.VariableNames{2} = 'NLac';              % adjust varnames
OUT2.Properties.VariableNames{1} = 'Lac';               % adjust varnames


%% STEP 3: check for outliers 3x std
for i = 1:length(OUT2.Lac(:,1)) 
    for ii = 3:length(OUT.Properties.VariableNames)
        ind = find(OUT.Lac ==i & (OUT{:,ii} > OUT2{i,ii}+3*STD(i,ii) | OUT{:,ii} < OUT2{i,ii} - 3*STD(i,ii)));
        OUT{ind,ii} = NaN;
    end
end

% recalculate herd averages, this time withoud outliers
OUT2 = unique(OUT.Lac); OUT2(:,2:length(OUT.Properties.VariableNames)) = 0;               % OUT2 will contain per lactation chars
for i = 1:max(OUT2(:,1))            
    OUT2(i,2) = sum(OUT.Lac == i);                      % number of lactations in ith lac for this herd
    OUT2(i,3:length(OUT.Properties.VariableNames)) = nanmean(OUT{OUT.Lac == i,3:length(OUT.Properties.VariableNames)},1);   % averages per lactation number
end
OUT2 = array2table(OUT2,'VariableNames',OUT.Properties.VariableNames); % array to table
OUT2.Properties.VariableNames{2} = 'NLac';              % adjust varnames
OUT2.Properties.VariableNames{1} = 'Lac';               % adjust varnames


%% STEP 3: also output herd-lactation corrected features

OUT3 = OUT;
for i = 1:length(OUT2.Lac(:,1))
    OUT3{OUT3.Lac==i,4:end} = OUT3{OUT3.Lac==i,4:end}./OUT2{i,4:end};
end
