function [cows,cowlac,SUM,DAY] = F1_RankResilience(DAY,varargin)
% This function calculates a ranking of dairy cows based on their
% estimated resilience, taking into account the fact they were (not)
% culled, reproduced and the time they got to get pregnant.
% 
% INPUTS:   DAY     'DailyData', containing following variables:
%                     - CowID       Unique Cow Identification number
%                                   --> ensure that a unique CowID is
%                                       provided for each of the cows in
%                                       the dataset!!
%                     - Lac         Lactation Number
%                     - Date        Date of the milking
%                     - Calving     Calving Date
%                     - DIM         Days in Milk
%                     - TMY         Total daily Milk Yield
%                     - BDay        BirthDay
%                     - (opt)SENSOR DATA: ECLF/ECLR/ECRF/ECRR/CowWeight/
%                                         MilkT/
%                     - (opt)DISEASE DATA
%                     - (opt)REPRODUCTION DATA (Ins/InsN)
%                     - (opt)ID     User defined Cow Identification Number
%                     - (opt)LacID  User defined Lactation ID Number
%
%           ACT     'Activity data'; this will be used to check whether 
%                   there is activity data and reproduction data available 
%                   for this lactation. If so, should contain following
%                   variables:
%                     - CowID       Unique Cow Identification number
%                                   --> ensure that a unique CowID is
%                                       provided for each of the cows in
%                                       the dataset!!
%                     - Lac         Lactation Number
%                     - Date        Date AND time of the measurement
%                     - Act         2 hourly activity data (might be
%                                   adjusted for other formats in a later 
%                                   stage)
%                     -(opt)Heat    Heat index
%                     -(opt)Att     Heat attention
%                     -(opt)Rum     Rumination value
%                     -(opt)HealthIndex HealthIndex
%
% OUTPUTS:  cows containing the overview of the unique cows and the
%                   rankings calculated for these (points)
%                   - CowID         Unique CowID
%                   - FirstLac      First Lactation number
%                   - LastLac       Last lactation number
%                   - UniqLac       Number of unique lactations
%                   - Start         StartDate for this cow
%                   - End           EndDate for this cow
%                   - Ndays         Number of days of data for this cow
%                   - Nmilk         Number of milkings available for this
%                   - avCI          Average CI
%                   - stdCI         Standard deviation of the CI
%                   - min/maxCI     Min and maximum Calving interval
%                   - RANKX         Points for ranking X%  
%                   - Incl1/Incl2   Included for the first analysis 4A or
%                                   for 4B
%
%           cowlac containing the data, scores, and individual information
%                   for each lactation for which data is available
%                   - CowID//FirstMilk/LastMilk/Start/End/Nmilk/NDays
%                   - Lac           Lactation number
%                   - Calving       Calving Date
%                   - NHE           Number of Health Events
%                   - NI            Number of inseminations
%                   - TMYd          Total milk deviation from lactation
%                                   average in %
%                   - Cull          Points for early (<100D) culling in
%                                   days: N days she was culled before 100 
%                                   DIM
%                   - CI            Calving interval for this lactation
%                   - CIdev         Calving interval deviation 
%                   - RANKX         Xth ranking (expressed in points)
%                   - LengthDryPeriod Length of the dry period according to
%                                   the data
%                   - Incl2         Including for second analysis
%
%           SUM     Summary of the herd with herein all information of the
%                   herd that was used to calculate the ranking on
%           SUM.CalvingIntervalHerd average/std/min/max
%           SUM.CalvingIntervalLac  average per lactation
%           SUM.LengthDryPeriod     average/std/min/max              
%           SUM.SensorData          table indication what and ho much
%                                   sensor data for each lactation is 
%                                   present
%           DAY                     input data, excluding measurments where
%           only one day is available or only NaN are present 
%
%
% STEP 1 Make overview of cows and cowlacs based on 'milkings' dataset and
%        calculate calving intervals
% STEP 2 Calculate rankings based on the following steps
%           RANK1:  LacN - CIdev 
%           RANK2:  LacN - CIdev - events/cull
%           RANK3:  LacN - CIdev -               305dY
%           RANK4:  LacN - CIdev - Nins
%           RANK5:  LacN - CIdev - events/cull - 305dY
%           RANK6:  LacN - CIdev - events/cull -          Nins
%           RANK7:  LacN - CIdev -               305dY  - Nins   
%           RANK8:  LacN - CIdev - events/cull - 305dY  - Nins
%           RANK9:         CIdev
%           RANK10:        CIdev - events (no cull)
% STEP 3 Summarize availability of sensor data for the cows and lactations
%        and the number of measurements and the availability of diseases.
%        Eg. when only mastitis data is present, only cows having udder
%        health problems will be disadvantaged
% STEP 4 Add exclusion criteria to the dataset - which cows are not further
%        considered in the ranking e.g. bc absence of sensor data or 
%        absence of enough measurements


% load testset to construct functions (delete in final function)
% % % clear variables
% % % close all
% % % clc
% % % load('C:\Users\u0084712\Box Sync\Documents\RAFT solutions\A1_Data\D09_data.mat');       % ACT1 DAY1 cows1 cowlac1
% % % DAY = DAY9; clear DAY9 cows9 cowlac9
% % % ACT = ACT9; clear ACT9


%% STEP 0: check whether you have activity data input
if nargin == 2
    ACT = varargin{1};
end


%% STEP 1: OVERVIEW of COWS and LACTATIONS, incl. CALVING INTERVALS
% Find column positions
VN1 = find(string(DAY.Properties.VariableNames) == 'CowID');    % col position of CowID
VN2 = find(string(DAY.Properties.VariableNames) == 'Lac');      % col position of Lac
VN3 = find(string(DAY.Properties.VariableNames) == 'Calving');  % col position of Calving
VN4 = find(string(DAY.Properties.VariableNames) == 'BDay');  % col position of Calving

% Calculate the average yield per lacation and per day
SUM.AverageMilkYield305 = ones(max(DAY.Lac),306);               % prepare array, row = lac, col = DIM
for i = 1:max(DAY.Lac)      % for all lactations in the dataset
    for j = 0:305           % for all days (0 to 305)
        ind = find(DAY.Lac == i & DAY.DIM == j & DAY.TMY ~= 0);
        SUM.AverageMilkYield305(i,j+1) = nanmean(DAY.TMY(ind));   % calculate mean for this DIM and lactation
        SUM.MedianMilkYield305(i,j+1) = nanmedian(DAY.TMY(ind));   % calculate mean for this DIM and lactation
   end
end


% plot average yield functions
% % % figure; plot(SUM.AverageMilkYield305','LineWidth',2)
% % % legend('1','2','3','4','5','6','7','8','9','10');
% % % xlim([0 306]);xlabel('DIM (days)');ylabel('Daily milk yield (kg)')
% % % title('Average daily yield for parity 1 to 10')

% Overview per lactation
cowlac = unique(DAY(:,[VN1 VN2 VN3]),'rows');   % Summary of unique rows [CowID Lac Calving]
cowlac = sortrows(cowlac,[ 1 3]);               % sort rows of cowlac summary
cowlac.BDate(:,1) = NaT;                        % Birth Date
cowlac.AFC(:,1) = NaN;
cowlac.FirstMilk(:,1) = NaN;                    % DIM first milking
cowlac.LastMilk(:,1) = NaN;                     % DIM last milking
cowlac.Start(:,1) = NaT;                        % Startdate for this lactation
cowlac.End(:,1) = NaT;                          % Enddate for this lactation
cowlac.Nmilk(:,1) = 0;                          % Number of milkings
cowlac.NDays(:,1) = 0;                          % Number of days data
cowlac.NHE(:,1) = 0;                            % Number of Health Events
cowlac.NI(:,1) = 0;                             % Number of inseminations
cowlac.TMYd(:,1) = 0;                           % Deviation of total average 305 yield
cowlac.Cull(:,1) = 0;                           % Number of days cow is culled before DIM 100
for i = 1:length(cowlac.Lac)
    try
        ind = find(DAY.Lac == cowlac.Lac(i) & strcmp(DAY.CowID(:,1),cowlac.CowID{i,1}) == 1 &  isnan(DAY.TMY) == 0);
    catch
        ind = find(DAY.Lac == cowlac.Lac(i) & DAY.CowID(:,1) == cowlac.CowID(i,1) &  isnan(DAY.TMY) == 0);
    end
    if isempty(ind) ==0 % if there's data available for this lactation, and TMY not 0 for all
        cowlac.FirstMilk(i,1) = min(DAY.DIM(ind));              % DIM first milking
        cowlac.LastMilk(i,1) = max(DAY.DIM(ind));               % DIM last milking
        cowlac.Start(i,1) = min(DAY.Date(ind));                 % Startdate for this lactation
        cowlac.End(i,1) = max(DAY.Date(ind));                   % Enddate for this lactation
        cowlac.Nmilk(i,1) = length(ind);                        % Number of milkings
        cowlac.NDays(i,1) = max(DAY.DIM(ind))-min(DAY.DIM(ind));% Number of days data
        if isempty(find(string(DAY.Properties.VariableNames) == 'Disease',1))==0    % check there's a column 'Disease' (optional)
            cowlac.NHE(i,1) = length((find(strcmp(DAY.Disease(ind,1),'') == 0)));   % number of health events detected
        end
        if isempty(find(string(DAY.Properties.VariableNames) == 'InsN',1))==0    % check there's a column 'Disease' (optional)
            cowlac.NI(i,1) = max([0 max(DAY.InsN(ind))]);                % Number of inseminations detected
        end
        % deviation of total yield, calculated for the days available for
        % this cow (max 305 days, but can also be 40-305 days e.g.)
        idx = find(DAY.DIM(ind) < 306);     % all indices for days smaller than 306
        cowlac.TMYd(i,1) = -(1-(nansum(DAY.TMY(ind(idx)))/nansum(SUM.AverageMilkYield305(cowlac.Lac(i),cowlac.FirstMilk(i)+1:min(cowlac.LastMilk(i)+1,306)))))*100;
        if cowlac.LastMilk(i) < 100 && cowlac.End(i) < max(DAY.Date) % make sure it is not just the end of the dataset 
            cowlac.Cull(i,1) = -(100-cowlac.LastMilk(i));            % penalties for LastMilk being before 100
        end
        cowlac.BDate(i,1) = DAY.BDay(ind(1));                   % add birth date
        if cowlac.Lac(i) == 1
            cowlac.AFC(i,1) = (datenum(cowlac.Calving(i,1))-datenum(cowlac.BDate(i,1)));
        end
    end
end
cowlac.TMYd(isnan(cowlac.TMYd) == 1,1) = 0;
clear idx i j ind

% find data for which startdate == enddate, and delete
% ind = find(datenum(cowlac.Start) == datenum(cowlac.End) | isnat(cowlac.Start)==1 | isnat(cowlac.End) == 1); % find all measurements
% cowlac(ind,:) = [];             % delete these measurements
% DAY = innerjoin(cowlac(:,[1 2 3]),DAY,'Keys',{'CowID','Lac','Calving'});
VN1 = find(string(DAY.Properties.VariableNames) == 'CowID');    % col position of CowID

% overview per cow
cows = unique(DAY(:,VN1));      % Summary of unique rows [CowID]
cows.AFC(:,1) = NaN;            % Age at first calving
cows.FirstLac(:,1) = NaN;       % First lactation for which data is available
cows.LastLac(:,1) = NaN;        % Last lactation for which data is available
cows.UniqLac(:,1) = NaN;        % Number of unique lactations for which data is available
cows.Start(:,1) = NaT;          % Start date of measurements for this cow
cows.End(:,1) = NaT;            % End date of measurments for this cow
cows.Ndays(:,1) = NaN;          % Number of days data is available for this cow
cows.Nmilk(:,1) = NaN;          % Number of milkings available for this cow
for i = 1:length(cows.CowID)
    try
        ind = find(strcmp(cowlac.CowID(:,1),cows.CowID{i,1}) == 1 & isnan(cowlac.FirstMilk)==0); % Overview based on cowlac, exclude lactations without data
    catch
        ind = find(cowlac.CowID(:,1) == cows.CowID(i,1)); % Overview based on cowlac, exclude lactations without data
    end
    
    if isempty(ind) == 0
        cows.FirstLac(i,1) = min(cowlac.Lac(ind));                   % First lactation
        cows.LastLac(i,1) = max(cowlac.Lac(ind));                    % Last lactation
        cows.UniqLac(i,1) = length(unique(cowlac.Lac(ind)));         % Number of unique lactations
        cows.Start(i,1) = min(cowlac.Start(ind));                    % First measurement of that cow
        cows.End(i,1) = max(cowlac.End(ind));                        % Last measurment of that cow
        cows.Ndays(i,1) =  max(datenum(cowlac.End(ind)))-min(datenum(cowlac.Start(ind))); % number of days included
        cows.Nmilk(i,1) = sum(cowlac.Nmilk(ind));                    % Number of milkings for this cow
        cows.AFC(i,1) = cowlac.AFC(ind(1));                          % Age at first calving
    end
end
clear i ind 

% Add calving interval to overview of lactations cowlac, add mean CI per
% cow to 'cows'
cowlac.CI(:,1) = NaN;                   % Calving Interval
cows.avCI(:,1) = NaN;                   % Average CI of this cow
cows.stdCI(:,1) = NaN;                  % Standard deviation of the CI of this cow
cows.minCI(:,1) = NaN;                  % minimum average CI of the cow
cows.maxCI(:,1) = NaN;                  % maximum average CI of the cow
for i = 1:length(cows.CowID)
    try
        ind = find(strcmp(cowlac.CowID(:,1),cows.CowID{i,1}) == 1); % find all lactations of cow with ID CowID(i)
    catch
        ind = find(cowlac.CowID(:,1) ==cows.CowID(i,1)); % find all lactations of cow with ID CowID(i)
    end
    if length(ind) == 1                 % of only one lactation is available for this cow
        cowlac.CI(ind) = 0;             % the CI is 0
        cows.avCI(i,1) = NaN;           % the average does not exist
        cows.stdCI(i,1) = NaN;          % the std does not exist
        cows.minCI(i,1) = NaN;          % the minCI does not exist
        cows.maxCI(i,1) = NaN;          % the maxCI does not exist
    else                                % if successive lactations are available for this cow
        starts = datenum(cowlac.Start(ind(1:length(ind)-1)) - cowlac.FirstMilk(ind(1:length(ind)-1)));  % start of the calvings
        newstarts = datenum(cowlac.Start(ind(2:end)) - cowlac.FirstMilk(ind(2:end)));                   % start of the next calvings
        CI = newstarts-starts;          % the calving intervals for this cow
        
        cowlac.CI(ind(1:end-1)) = CI;   % fill in the calving intervals
        cowlac.CI(ind(end)) = 0;        % the last calving interval is 0 (does not exist)
        cows.avCI(i,1) = mean(CI);      % average CI
        cows.stdCI(i,1) = std(CI);      % std CI within cow
        cows.minCI(i,1) = min(CI);      % min CI of this cow
        cows.maxCI(i,1) = max(CI);      % max CI of this cow
    end
end
clear starts newstarts cow CI ans S M L ind i cow

% add Herd mean CI
SUM.CalvingIntervalHerd = [nanmean(cows.avCI) nanmean(cows.stdCI) nanmean(cows.minCI) nanmean(cows.maxCI)]; % summary of the CI herd characteristics

% add Calving interval per lactation to calculate deviations
M = max(cowlac.Lac(cowlac.CI~=0));
for i = 1:M
    SUM.CalvingIntervalLac(i) = mean(cowlac.CI(cowlac.Lac == i & cowlac.CI~=0));
end

% add deviations to cowlac
cowlac.CIdev(cowlac.CI~=0,1) = -(cowlac.CI(cowlac.CI~=0)-SUM.CalvingIntervalLac(cowlac.Lac(cowlac.CI~=0,1))');       

% add average age at first calving to SUM
SUM.AFCHerd = [nanmean(cows.AFC) nanstd(cows.AFC) nanmin(cows.AFC) nanmax(cows.AFC)]; % summary of the CI herd characteristics


% add penalties for age at first calving, maximum 365 points subtracted
cowlac.AFCP(cowlac.AFC < 690,1) = 0.5*(cowlac.AFC(cowlac.AFC < 690,1) - 690);     % shorter than 23M
cowlac.AFCP(cowlac.AFC > 750,1) = max([(750 - cowlac.AFC(cowlac.AFC > 750,1)) -365*ones(length(cowlac.AFC(cowlac.AFC > 750,1)),1)],[],2);% longer than 25M
cowlac.AFCP(cowlac.AFC > 690 & cowlac.AFC < 750,1) = 750-cowlac.AFC(cowlac.AFC > 690 & cowlac.AFC < 750,1); % within 23 and 25 M

% add penalties for age at first calving based on the Herd average
cowlac.AFCS(:,1) = SUM.AFCHerd(1,1)-cowlac.AFC(:,1);

clear i M VN1 VN2 VN3


%% STEP 2: CALCULATE RANKING
% define weights
Wevents = -20;                                      % weight events
WlacN = 300;                                        % weight lacN -- 100 seems to be too low, 300 ensures higher lac have higher RS
WIns = [10 8 4 -2 -8:-6:-6*(max(cowlac.NI)-3)];     % weigth N inseminations
WAFC = 1;                                           % weigth for age at first calving

% add points based on deviations to cowlac
cowlac.RANK1(:,1) = NaN;    cows.RANK1(:,1) = NaN;     % RANK1:  LacN - CIdev - AFC 
cowlac.RANK2(:,1) = NaN;    cows.RANK2(:,1) = NaN;     % RANK2:  LacN - CIdev - AFC - events/cull
cowlac.RANK3(:,1) = NaN;    cows.RANK3(:,1) = NaN;     % RANK3:  LacN - CIdev - AFC -               305dY
cowlac.RANK4(:,1) = NaN;    cows.RANK4(:,1) = NaN;     % RANK4:  LacN - CIdev - AFC - Nins
cowlac.RANK5(:,1) = NaN;    cows.RANK5(:,1) = NaN;     % RANK5:  LacN - CIdev - AFC - events/cull - 305dY
cowlac.RANK6(:,1) = NaN;    cows.RANK6(:,1) = NaN;     % RANK6:  LacN - CIdev - AFC - events/cull -          Nins
cowlac.RANK7(:,1) = NaN;    cows.RANK7(:,1) = NaN;     % RANK7:  LacN - CIdev - AFC -               305dY  - Nins   
cowlac.RANK8(:,1) = NaN;    cows.RANK8(:,1) = NaN;     % RANK8:  LacN - CIdev - AFC - events/cull - 305dY  - Nins
cowlac.RANK9(:,1) = NaN;    cows.RANK9(:,1) = NaN;     % RANK9:         CIdev
cowlac.RANK10(:,1) = NaN;   cows.RANK10(:,1) = NaN;    % RANK10:        CIdev - events (no cull)

% The final all-inclusive score - RSCORE
cowlac.RSCORE(:,1) = NaN;   cows.RSCORE(:,1) = NaN;    % RESILIENCE SCORE BASED ON ALL AVAILABLE

VN4 = find(string(cows.Properties.VariableNames) == 'RANK1');       % col position of Rankings
VN5 = find(string(cowlac.Properties.VariableNames) == 'RANK1');     % col position of Rankings

cowlac.LengthDryPeriod(:,1) = NaN;                     % calculate length of dry period

for i = 1:length(cows.CowID)
    try
        ind = find(strcmp(cowlac.CowID(:,1),cows.CowID{i,1}) == 1);
    catch
        ind = find(cowlac.CowID(:,1) == cows.CowID(i,1));
    end
    cowlac.LengthDryPeriod(ind(1:end-1)) = datenum(cowlac.Calving(ind(2:end)))-datenum(cowlac.End(ind(1:end-1)));  % calculate max dry period for this farm
    
    for j = 1:length(ind)
        CumCIdev = sum(cowlac.CIdev(ind(1):ind(j)));     % cumulative CI deviation
        CumTMYd = sum(cowlac.TMYd(ind(1):ind(j)));       % cumulative TMY305 deviation
        CumEvent = sum(cowlac.NHE(ind(1):ind(j)));       % total number of events for this cow
        if cowlac.NI(ind(j)) == 0
            InsScore = 0;
        else
            InsScore = WIns(cowlac.NI(ind(j)));          % Insemination score 
        end
        AFCscore = nansum(cowlac.AFCP(ind));             % age at first calving score
        
        cowlac.RANK1(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev;   % RANK1:  LacN - CIdev 
        cowlac.RANK2(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + (Wevents*CumEvent + cowlac.Cull(ind(j)));  % RANK2:  LacN - CIdev - events/cull
        cowlac.RANK3(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + CumTMYd;     % RANK3:  LacN - CIdev - 305dY
        cowlac.RANK4(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + InsScore;        % RANK4:  LacN - CIdev - Nins
        cowlac.RANK5(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + (Wevents*CumEvent + cowlac.Cull(ind(j))) + CumTMYd;  % RANK5:  LacN - CIdev - events/cull - 305dY
        cowlac.RANK6(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + (Wevents*CumEvent + cowlac.Cull(ind(j))) + InsScore;        % RANK6:  LacN - CIdev - events/cull - Nins
        cowlac.RANK7(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + CumTMYd + InsScore;        % RANK7:  LacN - CIdev - 305dY       - Nins
        cowlac.RANK8(ind(j),1) = SUM.CalvingIntervalHerd(1) + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + (Wevents*CumEvent + cowlac.Cull(ind(j))) + InsScore;        % RANK8:  LacN - CIdev - events/cull - Nins   - 305dY
        cowlac.RANK9(ind(j),1) = SUM.CalvingIntervalHerd(1) + CumCIdev;                              % RANK9:  LacN - CIdev
        cowlac.RANK10(ind(j),1) = SUM.CalvingIntervalHerd(1) + CumCIdev + Wevents*CumEvent;          % RANK10: LacN - CIdev - events (no cull)
    
        cowlac.RSCORE(ind(j),1) = 500 + AFCscore + cowlac.Lac(ind(j),1)*WlacN + CumCIdev + (Wevents*CumEvent + cowlac.Cull(ind(j))) + InsScore + CumTMYd;
    end   
    
    % add final ranking to 'cows' overview - add 'standardized' number of
    % points for cows for which the data doesn't start from 1
    if cows.FirstLac(i) == 1
        cows(i,VN4:end) = cowlac(ind(j),VN5:end-1);           % add the final ranking, no standard points bc start = 1st lac
    else
        StPoints = SUM.CalvingIntervalHerd(1) + WlacN*(cows.FirstLac(i)-1);
        cows{i,VN4:end-2} = StPoints + cowlac{ind(j),VN5:end-3};  % add final ranking + standard points for lac without data
        cows{i,end-2:end} = SUM.CalvingIntervalHerd(1) + cowlac{ind(j),end-3:end-1};         % add RAFT ranking too
    end
end
clear ind i j StPoints VN1 VN2 VN3 VN4 VN5 CumCIdev CumEvent CumTMYd InsScore WIns WlacN Wevents ans


% calculate length of the dry period
SUM.LengthDryPeriod(1,1:4) = [nanmean(cowlac.LengthDryPeriod(cowlac.LengthDryPeriod>0)) nanstd(cowlac.LengthDryPeriod(cowlac.LengthDryPeriod>0)) min(cowlac.LengthDryPeriod(cowlac.LengthDryPeriod>0)) max(cowlac.LengthDryPeriod(cowlac.LengthDryPeriod>0))];                               % prepare finding the max length of the dry period for successive lactations


% plot rankings, sorted
% % % A = sortrows(cows, {'LastLac','RANK1'});        % sorted
% % % L = length(cows.CowID);             % 
% % % figure; plot(1:L,A.RANK1,1:L,A.RANK2,1:L,A.RANK3,1:L,A.RANK4,1:L,A.RANK5,1:L,A.RANK6,1:L,A.RANK7,1:L,A.RANK8,1:L,A.RANK9,1:L,A.RANK10,'LineWidth',2)
% % % legend({'1','2','3','4','5','6','7','8','9','10'},'Location','best');
% % % xlim([1 L])
% % % title('Points for ranking, sorted per lactation number and RANK1')
% % % xlabel('Cows in dataset')
% % % ylabel('Points for ranking')
% % % set(gca, 'FontSize',12)
% % % clear A L

%% STEP 3: ADD AVAILABILITY OF SENSOR DATA
% check whether there is sensor data (for the full lactation) available;
% the different possibilities implemented are:
%   1) Milk yield (but this is already in the summary
%   2) Electrical conductivity
%   3) SCC
%   4) Fat
%   5) Protein
%   6) Milk Temperature
%   7) Cow weights
%   8) Activity data (whatever format)

% % % % % L = length(cowlac.CowID);
% % % % % SUM.SensorData = cowlac(:,1:7);
% % % % % SUM.SensorData.MY(:,1:3) = ones(L,1)*[0 NaN NaN];    % yes, DIMstart, DIMend
% % % % % SUM.SensorData.EC(:,1:3) = ones(L,1)*[0 NaN NaN];    % yes, DIMstart, DIMend
% % % % % SUM.SensorData.SCC(:,1:3) = ones(L,1)*[0 NaN NaN];   % yes, DIMstart, DIMend
% % % % % SUM.SensorData.Fat(:,1:3) = ones(L,1)*[0 NaN NaN];   % yes, DIMstart, DIMend
% % % % % SUM.SensorData.Pro(:,1:3) = ones(L,1)*[0 NaN NaN];   % yes, DIMstart, DIMend
% % % % % SUM.SensorData.Tem(:,1:3) = ones(L,1)*[0 NaN NaN];   % yes, DIMstart, DIMend
% % % % % SUM.SensorData.Wei(:,1:3) = ones(L,1)*[0 NaN NaN];   % yes, DIMstart, DIMend
% % % % % SUM.SensorData.Act(:,1:3) = ones(L,1)*[0 NaN NaN];   % yes, DIMstart, DIMend
% % % % % for i = 1:length(SUM.SensorData.CowID)
% % % % %     ind = find(strcmp(DAY.CowID(:,1),SUM.SensorData.CowID{i,1}) == 1 & DAY.Lac == SUM.SensorData.Lac(i));  % find the data of this lactation in DAY
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'TMY', 1)) ==0                          % Only fill in if data available
% % % % %         SUM.SensorData.MY(i,1) = double(isempty(find(isnan(DAY.TMY(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.MY(i,1) ==1                                                               % Add start and end date
% % % % %             SUM.SensorData.MY(i,2) = DAY.DIM(ind(find(isnan(DAY.TMY(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.MY(i,3) = DAY.DIM(ind(find(isnan(DAY.TMY(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'ECLF', 1)) ==0                         % Only fill in if data available
% % % % %         SUM.SensorData.EC(i,1) = double(isempty(find(isnan(DAY.ECLF(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.EC(i,1) ==1                                                               % Add start and end date
% % % % %             SUM.SensorData.EC(i,2) = DAY.DIM(ind(find(isnan(DAY.ECLF(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.EC(i,3) = DAY.DIM(ind(find(isnan(DAY.ECLF(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'SCC', 1)) ==0                          % Only fill in if data available
% % % % %         SUM.SensorData.SCC(i,1) = double(isempty(find(isnan(DAY.SCC(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.SCC(i,1) == 1                                                             % Add start and end date
% % % % %             SUM.SensorData.SCC(i,2) = DAY.DIM(ind(find(isnan(DAY.SCC(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.SCC(i,3) = DAY.DIM(ind(find(isnan(DAY.SCC(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'Fat', 1)) ==0                          % Only fill in if data available
% % % % %         SUM.SensorData.Fat(i,1) = double(isempty(find(isnan(DAY.Fat(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.Fat(i,1) ==1                                                              % Add start and end date
% % % % %             SUM.SensorData.Fat(i,2) = DAY.DIM(ind(find(isnan(DAY.Fat(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.Fat(i,3) = DAY.DIM(ind(find(isnan(DAY.Fat(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'Protein', 1)) ==0                      % Only fill in if data available
% % % % %         SUM.SensorData.Pro(i,1) = double(isempty(find(isnan(DAY.Protein(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.Pro(i,1) ==1                                                              % Add start and end date
% % % % %             SUM.SensorData.Pro(i,2) = DAY.DIM(ind(find(isnan(DAY.Protein(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.Pro(i,3) = DAY.DIM(ind(find(isnan(DAY.Protein(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'MilkT', 1)) ==0                        % Only fill in if data available
% % % % %         SUM.SensorData.Tem(i,1) = double(isempty(find(isnan(DAY.MilkT(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.Tem(i,1) ==1                                                              % Add start and end date
% % % % %             SUM.SensorData.Tem(i,2) = DAY.DIM(ind(find(isnan(DAY.MilkT(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.Tem(i,3) = DAY.DIM(ind(find(isnan(DAY.MilkT(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if isempty(find(string(DAY.Properties.VariableNames) == 'CowWeight', 1)) ==0                    % Only fill in if data available
% % % % %         SUM.SensorData.Wei(i,1) = double(isempty(find(isnan(DAY.CowWeight(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.Wei(i,1) ==1                                                              % Add start and end date
% % % % %             SUM.SensorData.Wei(i,2) = DAY.DIM(ind(find(isnan(DAY.CowWeight(ind)) == 0,1,'first')),1);
% % % % %             SUM.SensorData.Wei(i,3) = DAY.DIM(ind(find(isnan(DAY.CowWeight(ind)) == 0,1,'last')),1);
% % % % %         end
% % % % %     end
% % % % %     if nargin == -2
% % % % %         ind = find(strcmp(ACT.CowID(:,1),SUM.SensorData.CowID{i,1}) == 1 & datenum(ACT.Date) > datenum(SUM.SensorData.Calving(i)) & datenum(ACT.Date) < datenum(SUM.SensorData.End(i))); 
% % % % %         SUM.SensorData.Act(i,1) = double(isempty(find(isnan(ACT.Act(ind)) == 0, 1))==0);
% % % % %         if SUM.SensorData.Act(i,1) ==1                                                              % Add start and end date
% % % % %             SUM.SensorData.Act(i,2) = floor(datenum(ACT.Date(ind(find(isnan(ACT.Act(ind)) == 0,1,'first')),1))-datenum(SUM.SensorData.Calving(i)));
% % % % %             SUM.SensorData.Act(i,3) = ceil(datenum(ACT.Date(ind(find(isnan(ACT.Act(ind)) == 0,1,'last')),1))-datenum(SUM.SensorData.Calving(i)));
% % % % %         end
% % % % %     end
% % % % % end
% % % % % clear i ind L M



%% STEP 4: INCLUSION and EXCLUSION CRITERIA

% for the further analysis, we need to decide which lactations can be
% included; based on the availability of sensor (at least milk yield) data
% for her whole productive life / at least two completed lactations
% Because the analysis will be split up in 2 different 'phases', I defined
% two different sets of selection criteria:
%   1) a strict 'completeness' criterion for the 'rankings' analysis
%       - includes only cows having their whole productive life available
%       in the dataset:
%           * starting from LACTATION 1 and DIM < 10
%           * ending before the end of the dataset (minus the dry period
%           correction of 120 days), thus date of last measurement is
%           before (end dataset - 120). 
%   2) a somewhat less strict completeness criterion for the longitudinal
%     'recalving' analysis
%       - includes all cows and lactations that have at least data 
%         available for 2 successive lactations, or evidence that she's
%         culled before end of dataset
% I've thus added two criteria: Incl1 is for the rankings analysis (only 
% cows for which the whole productive life is available = selection at cow 
% level) and Incl2 are for the 'recalving' analysis  (also cows not yet
% culled when dataset ends, but having enough data available to predict at 
% least once whether or not she will re-calve) = selection at cowlac level.

% % % % % cows.Incl1(:,1) = 1;    % data available for the cows whole productive life
% % % % % cows.Incl2(:,1) = 1;    % if (part of) the data is OK for the longitudinal re-calving analysis for this cow
% % % % % cowlac.Incl2(:,1) = 1;  % select at lactation level first!
% % % % % 
% % % % % for i = 1:length(cows.CowID)                                    % find all lactations of each cow
% % % % %     ind = find(strcmp(cowlac.CowID(:,1),cows.CowID{i,1}) == 1); 
% % % % %     Track = 0;
% % % % %     % Include 1 criterion completeness for ranking analysis
% % % % %     if cows.FirstLac(i) > 1
% % % % %         cows.Incl1(i) = 0;   % exclude if first lacs is not first parity
% % % % %     end
% % % % %     if cows.RANK5(i) < 0 
% % % % %         cows.Incl1(i) = 0;   % exclude if ranking is unreliable / error in lac
% % % % %         Track = 1;
% % % % %     end 
% % % % %     if cows.LastLac(i) ~= cows.UniqLac(i)
% % % % %         cows.Incl1(i) = 0;   % exclude if not all lactations contain data
% % % % %     end
% % % % %     if cowlac.FirstMilk(ind(1)) > 5
% % % % %         cows.Incl1(i) = 0;   % exclude if first milking of first lac later than DIM 5
% % % % %     end
% % % % %     if datenum(cowlac.End(ind(end))) > datenum(max(DAY.Date))-120
% % % % %         cows.Incl1(i) = 0;   % exclude if last milking is before end of dataset minus 'safety' hypothetical dry-off period of 120 days
% % % % %     end
% % % % %     
% % % % %     % Include 2 criterion at lactation level for 'longitudinal' re-calving analysis
% % % % %     for j = 1:length(ind)
% % % % %         if Track == 1
% % % % %             cowlac.Incl2(ind(j)) = 0;        % exclude if data is unreliable
% % % % %         elseif cowlac.FirstMilk(ind(j)) > 5  % if this is not the last lactation for this cow - include her! (at least one more calving date)
% % % % %             cowlac.Incl2(ind(j)) = 0;
% % % % %         elseif j ~= length(ind)
% % % % %             cowlac.Incl2(ind(j)) = 1;
% % % % %         elseif j == length(ind) && datenum(cowlac.End(ind(j))) > datenum(max(DAY.Date))-120  % if this is the last lactation for this cow & we are not sure she is culled before end of data - exclude her!   
% % % % %             cowlac.Incl2(ind(j)) = 0;
% % % % %         else                                % if she is culled before the end of the dataset - include her
% % % % %             cowlac.Incl2(ind(j)) = 1;
% % % % %         end        
% % % % %     end
% % % % % 
% % % % %     % Add second inclusion to cows.Incl2
% % % % %     cows.Incl2(i) = sum(cowlac.Incl2(ind)); %  if it needs to be 0/1, use: min([1 sum(cowlac.Incl2(ind))]);
% % % % % end
% % % % % clear idx i ind j ans A





