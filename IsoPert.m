function [OUT] = IsoPert(DIM,TMY,MOD,minloss,minlength,recpar,P,A)
% this function returns the different perturbations of a lactation, based
% on the DIM, TMY and the model given by MOD.
% All three inputs have the same size
% OUT will be an overview of the pertubarios containing following data:
%       1       number/order
%       2       start DIM
%       3       end DIM
%       4       maximum loss
%       5       
%
%
% if P == 1 -> absolute loss
% if P == 2 -> relative (%) loss
% STEP 1: calculate the residuals of TMY against MOD
% STEP 2: identify indices of negative residuals below minloss
% STEP 3: identify start and end of perturbation around minloss
% STEP 4: merge perturbations with overlapping end/start
% STEP 5: delete perturbations shorter than minlength
% STEP 6: calculate recovery - length (recpar = maxlength of recovery)
% STEP 7: outputs


% % % % MOD = Wood(p,DIM); 
% % % % minloss = -10;
% % % % minlength = 5;
% % % % recpar = 30;
% % % % OUT = [];
% % % % P = 2; A = 2;

% put robustfit warnings off
warning('off','stats:statrobustfit:IterationLimit')

% ensure negatives
minloss = -abs(minloss);


% STEP 1: calculatate the residuals
RES = TMY-MOD;                  % calculate residuals


% STEP 2: indices of perturbations below minloss
if P == 1                                                   % ABSOLUTE LOSS
    idx = find(RES < minloss);                              % minloss expressed as absolute values
elseif P == 2                                               % RELATIVE LOSS
    idx = find((RES./TMY)*100 < minloss);        % minloss expressed as absolute values
end

OUT(:,1) = TMY(idx);                                    % prepare output TMY
OUT(:,2) = DIM(idx);                                    % prepare output DIM


% STEP 3: start and end of each perturbation, DIM and index
for i = 1:length(idx)
    idx1 = find(RES > 0 & DIM < DIM(idx(i)), 1, 'last');% index of last above 0
    idx1 = idx1+1;                                      % index of first below 0
    if isempty(idx1) ==1; idx1 = 1; end                 % no index then is first
    OUT(i,3) = idx1;                                    % fill in index
    OUT(i,4) = DIM(idx1);                               % fill in DIM
    
    idx2 = find(RES > 0 & DIM > DIM(idx(i)), 1,'first');% index first above 0
    idx2 = min(idx2+1,length(DIM));                     % index of last below 0
    if isempty(idx2) ==1; idx2 = length(DIM); end       % no index then is last
    OUT(i,5) = idx2;                                    % fill in index
    OUT(i,6) = DIM(idx2);                               % fill in DIM
end

% STEP 4: merge perturbations with overlapping end/start
OUT(:,7) = 1;                                           % prepare OUT
for i = 1:length(idx)-1                                 % run through all perturbations
    if OUT(i+1,4) <= OUT(i,6)                           % if the previous overlaps
        OUT(i+1,7) = 0;                                 % then indicate to delete
    end
end
OUT(OUT(:,7)==0,:) = [];                                % delete

% STEP 5: calculate length of the perturbations & select
OUT(:,8) = OUT(:,6) - OUT(:,4);                         % calculate length
OUT(OUT(:,8) < minlength,:) = [];                       % delete short perturbations
    

if A == 1         % FIRST POSSIBILITY : TAKE INTO ACCOUNT ALL RESIDUALS OF PERTURBATION
    % STEP 6: calculate recovery parameter - length
    %               two-step procedure
    %               first step: calculate x from min to all
    %               second step: if intercept with x later than recpar:
    %               recalculate
    
    OUT(:,9:15) = NaN;                                  % prepare OUT
    for i = 1:length(OUT(:,1))                          % for all detected perturbations
        pertres = RES(DIM >= OUT(i,4) & DIM <= OUT(i,6));% the residuals of the Wood model within each perturbation
        smotres = smooth(pertres,length(pertres),'sgolay',3);% second order savitsky golay smoother
        [m,idx] = min(smotres);                         % min value of ith perturbation based on smoothed function
        m = min(m,RES(OUT(i,3)+idx(1)-1));              % chose between smoothed min and residual min to start regression from
        OUT(i,9) = m;                                   % fill in raw residual based on index of min smooth value
        
        OUT(i,10) = OUT(i,3)+idx(1)-1;                  % index of minimum
        OUT(i,11) = DIM(OUT(i,10));                     % DIM of minimum

%         m = min(RES(DIM >= OUT(i,4) & DIM <= OUT(i,6)));% min value of ith perturbation
%         OUT(i,9) = m;                                   % fill in min value
%         
%         ind = find(RES == m & DIM >= OUT(i,4) & DIM <= OUT(i,6), 1, 'first'); % index of minimum
%         OUT(i,10) = ind;                              % index of minimum
%         OUT(i,11) = DIM(ind);                         % DIM of index of minimum
  
    % Recovery phase
        PTc = RES(OUT(i,10):OUT(i,5))-OUT(i,9);         % perturbation - corrected
        DIMc = DIM(OUT(i,10):OUT(i,5))-DIM(OUT(i,10));  % DIM of the perturbation correct for min value
        TH = -OUT(i,9);                                 % threshold for crossing
        
        if length(PTc) > 3
            b = regress(PTc,DIMc);                      % slope of recovery after perturbation
            b = max(0.001,b);                           % positive slope required
            OUT(i,12) = min(recpar, TH/b);              % TH/b is moment that fitted line crosses '0' loss - number of days needed for recovery
            OUT(i,13) = TH*min(recpar, TH/b)/2;         % total 'losses' in this period - triangle surface of recovery period
        else
            OUT(i,12) = NaN;                            % No recovery data available
            OUT(i,13) = NaN;                            % Total 'losses' in this period = NaN
        end
        
    % Development phase
        PTc = RES(OUT(i,3):OUT(i,10))-m;                % perturbation - corrected
        DIMc = DIM(OUT(i,3):OUT(i,10))-DIM(OUT(i,10));  % DIM of the perturbation correct for min value
        TH = -m;                                        % Threshold for crossing
        
        if length(PTc) >= 2
            b = regress(PTc,DIMc);                      % Slope of recovery after perturbation
            b = min(-0.00001,b);                        % Negative slope required
            OUT(i,14) = min(OUT(i,11)-OUT(i,4), -TH/b); % TH/b is moment that fitted line crosses '0' loss - number of days needed for recovery
            OUT(i,15) = TH*OUT(i,14)/2;                 % Total 'losses' in this period - triangle surface of development period
        else
            OUT(i,14) = NaN;                            % Development is fast - 'unlimited'
            OUT(i,15) = m;                              % Total 'losses' in this period
        end
    end
    
    
elseif A == 2           % SECOND POSSIBILITY - NICS BROKEN STICK APPROACH
    
    for i = 1:length(OUT(:,1))                              % for all the perturbations
        pertres = RES(DIM >= OUT(i,4) & DIM <= OUT(i,6));   % the residuals of the Wood model within each perturbation
        if length(pertres)>3
            smotres = smooth(pertres,length(pertres),'sgolay',3);% third order savitsky golay smoother
        else
            smotres = smooth(pertres,length(pertres),'sgolay',2);% second
        end
        [m,idx] = min(smotres);                             % min value of ith perturbation based on smoothed function
        m = min(m,RES(OUT(i,3)+idx(1)-1));              % chose between smoothed min and residual min to start regression from

        OUT(i,9) = m;                                       % fill min
        
        OUT(i,10) = OUT(i,3)+idx(1)-1;                      % index of minimum
        OUT(i,11) = DIM(OUT(i,10));                         % DIM of index of minimum
                
        
        if OUT(i,10) < OUT(i,5)-2 && OUT(i,5)-OUT(i,10) >= 5% if DIM of the minimum is smaller than end of perturbation - 2 days    
            
        % Recovery phase    
            PTc = RES(OUT(i,10):OUT(i,5))-m;                % perturbation - corrected for value at minimum smooth
            DIMc = DIM(OUT(i,10):OUT(i,5))-DIM(OUT(i,10));  % DIM of the perturbation correct for min value
            TH = -m;                                        % threshold for cut regression line with pert
            
            tel = 1; SS = zeros(1,7);                       % prepare array
            for ii = 3:1:length(PTc)-3                      % idx of the window end of first line
                
                x1 = DIMc(1:ii);                            % x values of first regression line
                x2 = DIMc(ii+1:end);                        % x values of second regression line
                
                y1 = PTc(1:ii);                             % y values of second regression line
                y2 = PTc(ii+1:end);                         % y values of second regression line
                
                b1 = robustfit(x1,y1,'andrews',1.339,'off');% slope of regression line through first ii points // b1 = regress(y1,x1);   
                b1 = max(b1,0);                             % make sure it is positive
                b2 = robustfit(x2,y2,'andrews',1.339,'on'); % intercept and slope of regression line through second ii points //b2 = regress(y2,[ones(length(x2),1) x2]);
                
%                 IC = b2(1)/(b1(1)-b2(2));                   % intercept of the two lines
                
%                 ix1 = find(DIMc <  IC);                     % DIMc smaller than intercept of two lines
%                 ix2 = find(DIMc >= IC);                     % DIMc larger than intercept of two lines     
                ix1 = 1:ii;                                  % DIMc smaller than intercept of two lines
                ix2 = ii+1:length(DIMc);                     % DIMc larger than intercept of two lines     

                ss1 = sum((PTc(ix1) - b1(1)*DIMc(ix1)).^2);  % sums of squares for the first line
                ss2 = sum((PTc(ix2) - (b2(1)+b2(2)*DIMc(ix2))).^2); % sums of squared for the second line
                
                SS(tel,:) = [ii b1(1) b2(1) b2(2) ss1 ss2 ss1+ss2];% tel, slope-line1 intercept-line2 slope-line2 totalsumsofsquares
                
                tel = tel+1;                                % increase teller
            end
            indmin = find(SS(:,7) == min(SS(:,7)),1);       % find the minimum of total sums of squares
            b = SS(indmin,2);                               % slope of line producing minimum sums of squares
            
            OUT(i,12) = min(recpar, TH/b);                  % TH/b is moment that fitted line crosses '0' loss - number of days needed for recovery
            OUT(i,13) = TH*OUT(i,12)/2;                     % total 'losses' in this period - triangle surface of recovery period
            
        else
            OUT(i,12) = NaN;                                % No recovery data available
            OUT(i,13) = NaN;                                % total 'losses' in this period = NaN
        end 
        
        
% % % %                     figure; subplot(1,2,1); hold on;box on; title('Perturbation, 3rd order poly smoothing, regressions'); xlabel('DIM (days)'); ylabel('Residual milk yield (kg)')
% % % %                     subplot(1,2,1);
% % % %                     plot(DIM(OUT(i,3):OUT(i,5)), RES(OUT(i,3):OUT(i,5)),'o','LineWidth',2);
% % % %                     plot(DIM(OUT(i,3):OUT(i,5)), smotres,'LineWidth',2);
% % % %                     for ii = 1:5:length(SS(:,1))
% % % %                         subplot(1,2,1); hold on
% % % %                         plot(DIM(OUT(i,10):OUT(i,10)+SS(ii,1)),(OUT(i,9)*0:SS(ii,1))*SS(ii,2)+OUT(i,9),'b--')
% % % %                         plot(DIM(OUT(i,10)+SS(ii,1)+1:OUT(i,5)),((OUT(i,10)+SS(ii,1)+1:OUT(i,5))-(OUT(i,10)+SS(ii,1)+1))*SS(ii,4)+OUT(i,9)+(SS(ii,3)+SS(ii,4)*SS(ii,1)),'k--')
% % % %                     end
% % % %                     subplot(1,2,2); hold on;box on; title('Total sums of squares of bended regression'); xlabel('Points from breakpoint');ylabel('Sums of squares')
% % % %                     plot(SS(1:5:length(SS(:,1)),1),SS(1:5:length(SS(:,1)),7),'o','LineWidth',2)
% % % %                     plot(SS(:,1),SS(:,7),'-','LineWidth',2)

        
    % Development phase
        PTc = RES(OUT(i,3):OUT(i,10))-m;                % perturbation - corrected
        DIMc = DIM(OUT(i,3):OUT(i,10))-DIM(OUT(i,10));  % DIM of the perturbation correct for min value
        TH = -m;                                        % Threshold for crossing
        
        if length(PTc) >= 2
            b = regress(PTc,DIMc);                      % Slope of recovery after perturbation
            b = min(-0.00001,b);                        % Negative slope required
            OUT(i,14) = min(OUT(i,11)-OUT(i,4), -TH/b); % TH/b is moment that fitted line crosses '0' loss - number of days needed for recovery
            OUT(i,15) = TH*OUT(i,14)/2;                 % Total 'losses' in this period - triangle surface of development period
        else
            OUT(i,14) = NaN;                            % Development is fast - 'unlimited'
            OUT(i,15) = m;                              % Total 'losses' in this period
        end   
    end
end



% STEP 7: Prepare outputs
if isempty(OUT)==0
    OUT = [(1:length(OUT(:,1)))' OUT(:,1:6) OUT(:,8:15)];
    OUT = array2table(OUT,'VariableNames',{'Npert','TMY','DIM','IDXstart','DIMstart','IDXend','DIMend','Length','MinLos','IDXmin','DIMmin','RECO','Loss','DEVE','DLOSS'});
else
    OUT = array2table(zeros(1,15),'VariableNames',{'Npert','TMY','DIM','IDXstart','DIMstart','IDXend','DIMend','Length','MinLos','IDXmin','DIMmin','RECO','Loss','DEVE','DLOSS'});
end
    
% plot for examples
% figure;
% subplot(2,1,1); hold on; xlabel('DIM (days)'); ylabel('Milk yield (kg)'); box on;
% plot(DIM,TMY,'o-','LineWidth',2)
% plot(DIM,MOD,'-','LineWidth',2)
% if OUT.Npert(1) ~= 0
%     for i = 1:length(OUT.DIM)
%         plot(DIM(OUT.IDXstart(i):OUT.IDXend(i)), TMY(OUT.IDXstart(i):OUT.IDXend(i)),'ro-','LineWidth',2,'MarkerFaceColor','r','MarkerSize',3)
%     end
% end
% 
% subplot(2,1,2); hold on; xlabel('DIM (days)'); ylabel('Residual milk yield (kg)'); box on;
% plot(DIM,RES,'o-','LineWidth',2);
% plot(DIM,zeros(length(DIM),1),'-','LineWidth',2)
% if OUT.Npert(1) ~= 0
%     plot(OUT.DIMmin,OUT.MinLos,'ro','LineWidth',2,'MarkerFaceColor','r')
%     for i = 1:length(OUT.Loss)
%         plot([OUT.DIMmin(i) OUT.DIMmin(i)], [0 OUT.MinLos(i)],'m','LineWidth',2)
%         plot([OUT.DIMmin(i) OUT.DIMmin(i)+OUT.RECO(i)], [OUT.MinLos(i) 0],'m','LineWidth',2)
%         plot([OUT.DIMmin(i)-OUT.DEVE(i) OUT.DIMmin(i)], [0 OUT.MinLos(i)],'m','LineWidth',2)
%     end
% end


