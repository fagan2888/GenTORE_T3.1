function [CLSF, ALLCAS] = F3_ActivityFeatures(ACT,CLSF)
% This function calculates a number of sensor features for each cow and
% each cow lactation, based on the input 'data'. More specifically, the 
% function exists of X main parts, each agreeing with a certain the
% analysis for a certainsensor type; if there is time left, we will explore
% the possibility of also including a 'combined' set of sensors. 
% Ideally, our features are logically linked to resilience, and thus, 
% are ultimately useful to predict the resilience ranking. For this 
% function, I will start implementing sensor features calculated for
% 'completed' lactations; however in a second step we also want to have 
% a sort of monitoring tool to predict the cow starting a new lactation,
% which requires to include sensor features for part of the lactation, 
% more specifically, the beginning of a lactation (and the features of 
% the previous)
% The different types of sensor data used in this function are based on the
% 'LELY' datasets, and include:
%       MILK YIELD                  daily yields as a proxy for the cow's
%                                   whole physiological system
%       ACTIVITY                    as a proxy for fertility and general
%                                   health status
% INPUTS:   DATA                    dataset containing only data for the 
%                                   the already selected animals; i.e. for
%                                   the rankings-evaluation, this means
%                                   only cows with a complete history are
%                                   included in the analysis, selected on
%                                   the output variable 'Incl' from F1
%                                   at least following variables included:
%                                       CowID, Lac, Calving, Date, DIM
%           varargin             OR activity data table, containing
%                                   variables included:
%                                       CowID, date, Act
%                                   optional:
%                                       Heat, Att, Rum
%                                OR sensor feature indication (array FEAT)
%                                       1 = MILK YIELD
%                                       2 = ACTIVITY
%                                       3 = WEIGHT
%                                       4 = ELECTRICAL CONDUCTIVITY
%                                       5 = SCC
%                                       6 = FAT/PROTEIN
%                                       7 = TEMPERATURE
%                                OR both activity data table AND feature
%                                   indication, with feat indication = last
%                                   argument
%
% STEP 1: check which features need to be calculated, either based on the
%         availablity of data, or based on the array FEAT and give warning
%         or error when the data available does not match the features to
%         be calculated
% STEP 2: define CSF = 'Cows Sensor Features' and 
%                CLSF = 'Cow Lactation Sensor Features' 
%         - depending on FEAT, features are added to these output variables
% STEP 3: calculate MILK YIELD features
%               -
%               -
% STEP 4: calculate ACTIVITY features 
%               -
%               -


%% STEP 0: load data to construct function
% % % % load('D1_trainingset') % act, cowlac, cows, data, Overview, SUM
% % % % clear cows cowlac Overview SUM  % this are the previous rankings
% % % % % load new rankings
% % % % load('C:\Users\u0084712\Box Sync\Documents\RAFT solutions\A4_results\D_rankings2.mat')
% % % % 
% ACT = ACT_Pyrland; CLSF = CLSF1(CLSF1.FarmID == 21,[1 2 3 4 8 9]);

ALLCAS=0;
%% STEP 1: check data vs feature calculation
% only if the array FEAT is given as an input, this is important
% check number of input arguments > 0



%% STEP 2: Define CSF and CLSF the features should be added to these tables
% find needed variables


%% STEP 4: SENSOR FEATURES -- ACTIVITY
% if ACTIVITY is included as sensor data, this part of the function will
% calculate sensor features for activity characteristics
% I would like to include both 'black box' features and 'non-black box'
% features, the latter for which we think of them being significant for
% resilience. For example, regularity of the activity peaks, timing of the
% first peak and peak height might be indicative over a cows' reproduction
% performance and health. However, also other features are possible.
% We will start the calculation of the features based on daily activity
% aggregates, while we will also consider the variability of the activity
% data. 
% The activity sensor data needs two additional preprocessing steps before
% we can start working on it/
%       1) add DAYS in MILK
%       2) add Lactation Number

% predefine the different cols for filling in feature characteristics
% CAT1 - general activity data charactistics
CLSF.ACStart(:,1) = NaN;            % Start DIM of the activity data
CLSF.ACEnd(:,1) = NaN;              % End DIM of the activity data
CLSF.ACDailyMean(:,1) = NaN;        % Mean total activity for that cow
CLSF.ACNumberZero(:,1) = NaN;       % Number of days the activity sum is zero
CLSF.ACPercZero(:,1) = NaN;         % Percentage of zero activity 
CLSF.ACSkewness(:,1) = NaN;         % Skewness of the activity data
CLSF.ACOverallVariance(:,1) = NaN;  % Average within-day variance in activity for that cow
CLSF.ACLacSlope(:,1) =  NaN;        % Slope of the activity

% CAT2 - identification of peaks & smoothing characteristics
% LEVEL 1 = raw data compared to median smoothed 4 day window
% this level should represent the 'extreme' peaks compared to the
% normal averages, lasting 1 day, and thus is what we should normally
% discriminate as oestrus-related features
CLSF.ACL1Npeaks(:,1) = NaN;         % Number of high peaks in the raw data
CLSF.ACL1PeakHeiAv(:,1) = NaN;      % Average peak height
CLSF.ACL1PeakHeiMax(:,1) = NaN;     % Maximum peak height
CLSF.ACL1Npeaks150D(:,1) = NaN;     % Number of high peaks in the first 150 days of data
CLSF.ACL1AvPeakInter(:,1) = NaN;    % Average interval between successive peaks
CLSF.ACL1StdPeakInter(:,1) = NaN;   % Standard deviation on the interval between successive peaks
CLSF.ACL1AvPeakInter150D(:,1) = NaN;% Average peak interval in the first 150 days of data
CLSF.ACL1StdPeakInter150D(:,1) = NaN;% Standard deviation of the between-peak interval in the first 150 days of data
CLSF.ACL1FirstPeak(:,1) = NaN;      % DIM of the first day

% LEVEL 2 = median smoothed 4D vs median smoothed 20D
% this level represents the variability in average level of activity,
% and might thus be related to the temperament and health events of the
% cows. Both peaks and valleys are detected, and we also take
% regularity into account
CLSF.ACL2NPosPeaks(:,1) = NaN;      % Number of positive peaks (more than average activity)
CLSF.ACL2PosPeakHeiAv(:,1) = NaN;   % Average height of the peaks
CLSF.ACL2PosPeakHeiMax(:,1) = NaN;  % Maximum height of the peaks
CLSF.ACL2PosAvPeakInter(:,1) = NaN; % Average interval between peaks
CLSF.ACL2PosStdPeakInter(:,1) = NaN;% Standard deviation on the interval between peaks
CLSF.ACL2PosAvPeakWidth(:,1) = NaN; % Width of the peaks (how long there is more or less activity than 'normal')
CLSF.ACL2NPosPeaks150D(:,1) = NaN;  % Number of positive peaks (more than average activity) in first 150 days
CLSF.ACL2PosPeaksInt150D(:,1) = NaN;% Average interval of the peaks in the first 150 D
CLSF.ACL2PosPeaksIntStd150D(:,1) = NaN;% Standard deviation on interval of the peaks in the first 150 D

CLSF.ACL2NNegPeaks(:,1) = NaN;      % Number of valleys (less than average activity)
CLSF.ACL2NegPeakHeiAv(:,1) = NaN;   % Average valley depth
CLSF.ACL2NegPeakHeiMax(:,1) = NaN;  % Maximum valley depth
CLSF.ACL2NegAvPeakInter(:,1) = NaN; % Average interval between valley
CLSF.ACL2NegStdPeakInter(:,1) = NaN;% Standard deviation on the interval between valleys
CLSF.ACL2NegAvPeakWidth(:,1) = NaN; % Width of the valleys (how long there is more or less activity than 'normal')
CLSF.ACL2NNegPeaks150D(:,1) = NaN;  % Number of valleys (more than average activity) in first 150 days
CLSF.ACL2NegPeaksInt150D(:,1) = NaN;% Average interval of the valleys in the first 150 D
CLSF.ACL2NegPeaksIntStd150D(:,1) = NaN;% Standard deviation on interval of the valleys in the first 150 D

CLSF.ACL2NAllPeaks(:,1) = NaN;      % Total number of bumps in the data
CLSF.ACL2AllAvPeakInter(:,1) = NaN; % Average peak interval
CLSF.ACL2AllStdPeakInter(:,1) = NaN;% Standard deviation on the peak intervals
CLSF.ACL2AllAvPeakWidth(:,1) = NaN; % Average bumps widths 
CLSF.ACL2NAllPeaks150D(:,1) = NaN;  % Number of bumps in first 150 days
CLSF.ACL2AllPeaksInt150D(:,1) = NaN;% Average interval of the bumps in the first 150 D
CLSF.ACL2AllPeaksIntStd150D(:,1) = NaN;% Standard deviation on interval of the bumps in the first 150 D


% %     % LEVEL 3 = median smoothed 20D vs intercept slope linear model
% %     % This features show the general behavior of the activity time series,
% %     % and represent the broader changes over time. 
% %     CLSF.ACL3NPosPeaks(:,1) = NaN;      % Number of periods more than normal activity 
% %     CLSF.ACL3PosPeakHeiAv(:,1) = NaN;   % Height of the deviation
% %     CLSF.ACL3PosPeakHeiMax(:,1) = NaN;  % maxim
% %     CLSF.ACL3PosAvPeakInter(:,1) = NaN;
% %     CLSF.ACL3PosStdPeakInter(:,1) = NaN;
% %     CLSF.ACL3PosAvPeakWidth(:,1) = NaN;
% %     
% %     CLSF.ACL3NNegPeaks(:,1) = NaN;
% %     CLSF.ACL3NegPeakHeiAv(:,1) = NaN;
% %     CLSF.ACL3NegPeakHeiMax(:,1) = NaN;
% %     CLSF.ACL3NegAvPeakInter(:,1) = NaN;
% %     CLSF.ACL3NegStdPeakInter(:,1) = NaN;
% %     CLSF.ACL3NegAvPeakWidth(:,1) = NaN;
% %     
% %     CLSF.ACL3NAllPeaks(:,1) = NaN;
% %     CLSF.ACL3AllAvPeakInter(:,1) = NaN;
% %     CLSF.ACL3AllStdPeakInter(:,1) = NaN;
% %     CLSF.ACL3AllAvPeakWidth(:,1) = NaN;


% predefine length of ALL CAS containing the daily activity data
ALLCAS = zeros(floor(length(ACT.CowID)/12),7);     
tel = 0;        % teller to fill in ALLCAS with summarized activity data

if sum(contains(CLSF.Properties.VariableNames,'CowID'))==0
    CLSF.CowID = CLSF.CowId;
end
    
for i = 1:length(CLSF.CowID)    % for all the cows, find and summarize the activity data
    try
        ind = find(strcmp(ACT.CowID(:,1),CLSF.CowID{i,1}) == 1 & datenum(ACT.Date(:,1))>= datenum(CLSF.Start(i)) & datenum(ACT.Date(:,1))<= datenum(CLSF.End(i)));    % find all data of cow and lac i
    catch
        ind = find(ACT.CowID(:,1)==CLSF.CowID(i,1) & datenum(ACT.Date(:,1))>= datenum(CLSF.Startdate(i)) & datenum(ACT.Date(:,1))<= datenum(CLSF.Enddate(i)));    % find all data of cow and lac i
    end
    % first step = set up new dataset which 'summarizes' the different
    % activity measures from the raw data
    if isempty(ind) == 0 && length(ind) > 24*10   % check whether we have activity data available for this cow lactation
        CowAct = ACT(ind,:);            % all activity data of this cow
        CowAct.Act = CowAct.Activity;
        CowAct.DIM = datenum(CowAct.Date)-datenum(CLSF.Startdate(i))+0.0001;     % add the days in milk to the data
        CowAct.Lac(:,1) = CLSF.Lac(i);

        CAS = [];           % cow activity summary
        CAS(:,1) = accumarray(ceil(CowAct.DIM),CowAct.Act,[],@nansum);     % we decided to work with raw 2-hourly measures,not with the attentions/rumination
        CAS(:,2) = accumarray(ceil(CowAct.DIM),CowAct.Act,[],@nanmean);    % accumulation of the activity data per day, using sum/mean/... 
        CAS(:,3) = accumarray(ceil(CowAct.DIM),CowAct.Act,[],@nanstd);
        CAS(:,4) = accumarray(ceil(CowAct.DIM),CowAct.Act,[],@nanvar);
        CAS(:,5) = accumarray(ceil(CowAct.DIM),CowAct.Act,[],@nanmin);
        CAS(:,6) = accumarray(ceil(CowAct.DIM),CowAct.Act,[],@nanmax);
        CAS = array2table(CAS, 'VariableNames',{'Actsum','Actmean','Actstd','ActVar','Actmin','Actmax'});
        CAS.DIM(:,1) = (1:length(CAS{:,1}))';


%             A1 = figure('Units','centimeters','OuterPosition',[1 1 25 18]); subplot(2,1,1); hold on; box on; xlabel('DIM (days)');ylabel('Daily sum activity')
%             plot(CAS.DIM,CAS.Actsum,'o-','MarkerSize',4,'MarkerFaceColor',[0.39 0.47 0.64],'LineWidth',1.5,'Color',[0.39 0.47 0.64])
%             xlim([0 max(CAS.DIM)])

        b = regress(CAS.Actsum,[ones(length(CAS.DIM),1),CAS.DIM]);
%                         plot(CAS.DIM,b(1)+b(2)*CAS.DIM,'b-','LineWidth',2)


        CLSF.ACStart(i,1) = CAS.DIM(find(CAS.Actsum~=0,1));         % Start DIM of the activity data
        CLSF.ACEnd(i,1) = CAS.DIM(find(CAS.Actsum~=0,1,'last'));    % End DIM of the activity data
        CLSF.ACDailyMean(i,1) = nanmean(CAS.Actsum);                % Mean total activity for that cow
%                         plot(CAS.DIM,ones(length(CAS.Actsum),1)*nanmean(CAS.Actsum),'--','LineWidth',2,'Color', [0 0.75 0.75])

        CLSF.ACNumberZero(i,1) = length(find(CAS.Actsum ==0));      % Number of days the activity sum is zero
        CLSF.ACPercZero(i,1) = length(find(CAS.Actsum==0))/length(CAS.DIM);      % Percentage of non-zero activity 
        CLSF.ACSkewness(i,1) = skewness(CAS.Actsum);                % Skewness of the activity data
        CLSF.ACOverallVariance(i,1) = nanmean(CAS.ActVar);          % Average within-day variance in activity for that cow
        CLSF.ACLacSlope(i,1) = b(2);        % Slope of the activity

        % CAT2 - identification of peaks

        % first step is to find a way to identify the peaks, we need to
        % compare this agains the average of the surrounding days, so
        % idealiter we calculate something as a 'leave one out' moving
        % average
        yy1 = movmedian(CAS.Actsum,4);
        Ryy1 = CAS.Actsum - yy1; 
        ma1 = movmad(Ryy1,7);
        Ryy1(:,2) = 4*ma1(:,1);
        Ryy1(:,3) = 0;
        Ryy1(Ryy1(:,1)> Ryy1(:,2) & Ryy1(:,1) > 0.4*max(Ryy1(:,1)),3) = 1; 
        Ryy1(:,4) = CAS.DIM;

        yy2 = movmedian(CAS.Actsum,20);

        Ryy2 = yy1 - yy2;
        a1 = []; a2 = [];
        if length(Ryy2) > 5
            [a1(:,1), a1(:,2), a1(:,3), a1(:,4)]= findpeaks(Ryy2,CAS.DIM,'MinPeakProminence',0.2*max(Ryy2),'MinPeakHeight',0.2*max(Ryy2),'MinPeakDistance',min(5,length(Ryy2(:,1))-2));%find(Ryy1(:,1) > Ryy1(:,2) & Ryy1(:,1) > 200);
            [a2(:,1), a2(:,2), a2(:,3), a2(:,4)]= findpeaks(-Ryy2,CAS.DIM,'MinPeakProminence',0.2*max(-Ryy2),'MinPeakHeight',0.2*max(-Ryy2),'MinPeakDistance',min(5,length(Ryy2(:,1))-2));%find(Ryy1(:,1) > Ryy1(:,2) & Ryy1(:,1) > 200);
        end
% %             Ryy3 = yy2 - (b(1)+b(2)*CAS.DIM); b1 = []; b2 = [];
% %             [b1(:,1), b1(:,2), b1(:,3), b1(:,4)]= findpeaks(Ryy3,CAS.DIM,'MinPeakProminence',30,'MinPeakDistance',10);
% %             [b2(:,1), b2(:,2), b2(:,3), b2(:,4)]= findpeaks(-Ryy3,CAS.DIM,'MinPeakProminence',30,'MinPeakDistance',10);


% % % % %              F5 = figure('Units','centimeters','OuterPosition',[1 1 26 18]); subplot(3,1,1); hold on; box on; xlabel('DIM (days)'); ylabel('Median corrected activity'),  xlim([0 max(CAS.DIM)])
% % % % %                         hold on; 
% % % % %                         if CLSF.ACNumberZero(i,1)>0; plot(CAS.DIM(CAS.Actsum == 0),CAS.Actsum(CAS.Actsum==0),'ro','MarkerSize',6,'MarkerFaceColor','r'); end
% % % % %                         plot(CAS.DIM,CAS.Actsum,'o-','MarkerSize',4,'MarkerFaceColor',[0.39 0.47 0.64],'LineWidth',1.5,'Color',[0.39 0.47 0.64])
% % % % %                         plot(CAS.DIM, yy1,'LineWidth',2)
% % % % %                         plot(CAS.DIM, yy2,'LineWidth',2)
% % % % %             
% % % % %                         subplot(3,1,2); hold on; box on; xlabel('DIM (days)'); ylabel('Median corrected activity'),  xlim([0 max(CAS.DIM)])
% % % % %                         plot(CAS.DIM, Ryy1(:,1),'o-','MarkerSize',4,'MarkerFaceColor',[0.6 0.2 0],'LineWidth',1.5,'Color',[0.6 0.2 0])
% % % % %                         plot(CAS.DIM(Ryy1(:,3) == 1),Ryy1(Ryy1(:,3) == 1,1),'rs','MarkerFaceColor','r')
% % % % %                         
% % % % %                         subplot(3,1,3); hold on; box on; xlabel('DIM (days)'); ylabel('Residual 20-4 days median smoothed activity'), xlim([0 max(CAS.DIM)])
% % % % %                         plot(CAS.DIM, Ryy2(:,1),'o-','MarkerSize',4,'MarkerFaceColor',[0.93 0.69 0.13],'LineWidth',1.5,'Color',[0.93 0.69 0.13])
% % % % %                         plot(a1(:,2),a1(:,1),'rs','MarkerFaceColor','r')
% % % % %                         plot(a2(:,2),-a2(:,1),'rs','MarkerFaceColor','r')

% %                         subplot(4,1,4); hold on; box on; xlabel('DIM (days)'); ylabel('Residual 20D median smoothed - lin trend activity'), xlim([0 max(CAS.DIM)])
% %                         plot(CAS.DIM, Ryy3(:,1),'o-','MarkerSize',4,'MarkerFaceColor','b','LineWidth',1.5,'Color','b')
% %                         plot(b1(:,2),b1(:,1),'rs','MarkerFaceColor','r')
% %                         plot(b2(:,2),-b2(:,1),'rs','MarkerFaceColor','r')


        % Now, characterize the peaks at each level in the data
        % LEVEL 1 = raw data compared to median smoothed 4 day window
        if isempty(Ryy1(:,1)) == 0
            CLSF.ACL1Npeaks(i,1) = sum(Ryy1(:,3));
            if CLSF.ACL1Npeaks(i,1) ~= 0
                CLSF.ACL1PeakHeiAv(i,1) = mean(Ryy1(Ryy1(:,3)==1,1));
                CLSF.ACL1PeakHeiMax(i,1) = max(Ryy1(Ryy1(:,3)==1,1));
                CLSF.ACL1Npeaks150D(i,1) = sum(Ryy1(Ryy1(:,4)<150,3));
                CLSF.ACL1AvPeakInter(i,1) = mean(diff(Ryy1(Ryy1(:,3)==1,4)));
                CLSF.ACL1StdPeakInter(i,1) = std(diff(Ryy1(Ryy1(:,3)==1,4)));
                CLSF.ACL1AvPeakInter150D(i,1) = mean(diff(Ryy1(Ryy1(:,3)==1&Ryy1(:,4)<150,4)));
                CLSF.ACL1StdPeakInter150D(i,1) = std(diff(Ryy1(Ryy1(:,3)==1&Ryy1(:,4)<150,4)));
                CLSF.ACL1FirstPeak(i,1) = Ryy1(find(Ryy1(:,3)==1,1),4);
            else
                CLSF.ACL1PeakHeiAv(i,1) = 10;
                CLSF.ACL1PeakHeiMax(i,1) = 10;
                CLSF.ACL1Npeaks150D(i,1) = 0;
                CLSF.ACL1AvPeakInter(i,1) = 0;
                CLSF.ACL1StdPeakInter(i,1) = 0;
                CLSF.ACL1AvPeakInter150D(i,1) = 0;
                CLSF.ACL1StdPeakInter150D(i,1) = 0;
                CLSF.ACL1FirstPeak(i,1) = 300;
            end
        end
        % LEVEL 2 = median smoothed 4D vs median smoothed 20D
        if isempty(a1)==0
            CLSF.ACL2NPosPeaks(i,1) = length(a1(:,1));
            CLSF.ACL2PosPeakHeiAv(i,1) = mean(a1(:,4));
            CLSF.ACL2PosPeakHeiMax(i,1) = max(a1(:,4));
            CLSF.ACL2PosAvPeakInter(i,1) = mean(diff(a1(:,2)));
            CLSF.ACL2PosStdPeakInter(i,1) = std(diff(a1(:,2)));
            CLSF.ACL2PosAvPeakWidth(i,1) = mean(a1(:,3));
            CLSF.ACL2NPosPeaks150D(i,1) = length(a1(a1(:,2)<150,1));  
            CLSF.ACL2PosPeaksInt150D(i,1) = mean(diff(a1(a1(:,2)<150,2)));
            CLSF.ACL2PosPeaksIntStd150D(i,1) = std(diff(a1(a1(:,2)<150,2)));
        end

        if isempty(a2) == 0
            CLSF.ACL2NNegPeaks(i,1) = length(a2(:,1));
            CLSF.ACL2NegPeakHeiAv(i,1) = mean(a2(:,4));
            CLSF.ACL2NegPeakHeiMax(i,1) = max(a2(:,4));
            CLSF.ACL2NegAvPeakInter(i,1) = mean(diff(a2(:,2)));
            CLSF.ACL2NegStdPeakInter(i,1) = std(diff(a2(:,2)));
            CLSF.ACL2NegAvPeakWidth(i,1) = mean(a2(:,3));
            CLSF.ACL2NNegPeaks150D(i,1) = length(a2(a2(:,2)<150,1));   
            CLSF.ACL2NegPeaksInt150D(i,1) = mean(diff(a2(a2(:,2)<150,2)));
            CLSF.ACL2NegPeaksIntStd150D(i,1) = std(diff(a2(a2(:,2)<150,2)));
        end

        if isempty(a2)==0 && isempty(a1) ==0
            a3 = sortrows([a1;a2],2);
            CLSF.ACL2NAllPeaks(i,1) = length(a3(:,1));
            CLSF.ACL2AllAvPeakInter(i,1) = mean(diff(a3(:,2)));
            CLSF.ACL2AllStdPeakInter(i,1) = std(diff(a3(:,2)));
            CLSF.ACL2AllAvPeakWidth(i,1) = mean(a3(:,3));
            CLSF.ACL2NAllPeaks150D(i,1) = length(a3(a3(:,2)<150,1));  
            CLSF.ACL2AllPeaksInt150D(i,1) = mean(diff(a3(a3(:,2)<150,2)));
            CLSF.ACL2AllPeaksIntStd150D(i,1) = std(diff(a3(a3(:,2)<150,2)));
        end





% %             % LEVEL 3 = median smoothed 20D vs intercept slope linear model
% %             if isempty(b1)==0
% %                 CLSF.ACL3NPosPeaks(i,1) = length(b1(:,1));
% %                 CLSF.ACL3PosPeakHeiAv(i,1) = mean(b1(:,4));
% %                 CLSF.ACL3PosPeakHeiMax(i,1) = max(b1(:,4));
% %                 CLSF.ACL3PosAvPeakInter(i,1) = mean(diff(b1(:,2)));
% %                 CLSF.ACL3PosStdPeakInter(i,1) = std(diff(b1(:,2)));
% %                 CLSF.ACL3PosAvPeakWidth(i,1) = mean(b1(:,3));
% %             end
% %             if isempty(b2) == 0
% %                 CLSF.ACL3NNegPeaks(i,1) = length(b2(:,1));
% %                 CLSF.ACL3NegPeakHeiAv(i,1) = mean(b2(:,4));
% %                 CLSF.ACL3NegPeakHeiMax(i,1) = max(b2(:,4));
% %                 CLSF.ACL3NegAvPeakInter(i,1) = mean(diff(b2(:,2)));
% %                 CLSF.ACL3NegStdPeakInter(i,1) = std(diff(b2(:,2)));
% %                 CLSF.ACL3NegAvPeakWidth(i,1) = mean(b2(:,3));
% %             end
% %             if isempty(b2)==0 && isempty(b1) ==0
% %                 b3 = sortrows([b1;b2],2);
% %                 CLSF.ACL3NAllPeaks(i,1) = length(b3(:,1));
% %                 CLSF.ACL3AllAvPeakInter(i,1) = mean(diff(b3(:,2)));
% %                 CLSF.ACL3AllStdPeakInter(i,1) = std(diff(b3(:,2)));
% %                 CLSF.ACL3AllAvPeakWidth(i,1) = mean(b3(:,3));
% %             end




        % I'd like to add still something about act in the non-preg
        % period
%             CLSF.ACNHighVarDay(i,1) = NaN;        % Number of high variance days

        ALLCAS(tel+1:tel+length(CAS.DIM),:) = CAS{:,:};  % fill in ALLCAS
        tel = tel+length(CAS.DIM);                  % increase teller

% % % % %             saveas(F5,['path\'  'ActFeat_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.jpg'],'jpg')
% % % % %             saveas(F5,['path\'  'MActFeat_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.fig'],'fig')

    end
end


ALLCAS = array2table(ALLCAS,'VariableNames',{'Actsum','Actmean','Actstd','Actvar','Actmin','Actmax','DIM'});
ALLCAS(ALLCAS.DIM == 0 & ALLCAS.Actsum == 0 & ALLCAS.Actstd == 0 & ALLCAS.Actmax == 0,:) = [];      % delete empty rows
