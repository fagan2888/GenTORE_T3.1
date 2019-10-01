function [CSF, CLSF, ALLCAS] = F2_SensorFeatures(DATA,varargin)
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
%       ELECTRICAL CONDUCTIVIY      four quarters separately, might need
%                                   extension to 1 combined measure as a
%                                   proxy for udder health
%       BODY WEIGHTS                as a proxy for feed intake and
%                                   production (E in/E out)
%       TEMPERATURE                 as a proxy for infection status (?)
%       FAT/PROT                    as a proxy for energy balance (?)
%       SCC                         as a proxy for udder health(?)
%
% INPUTS:   data                    dataset containing only data for the 
%                                   the already selected animals; i.e. for
%                                   the rankings-evaluation, this means
%                                   only cows with a complete history are
%                                   included in the analysis, selected on
%                                   the output variable 'Incl' from F1
%                                   at least following variables included:
%                                       CowID, Lac, Calving, Date, DIM
%                                   optional:
%                                       TMY, EC, CowWeight, MilkT, SCC, 
%                                       Fat, Protein
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
% % % % DATA = data.dataset_17;   % select the data of dataset 1 to work with
% % % % ACT = act.dataset_17;     % select activity data of dataset 1 to work with 
% % % % 

ALLCAS=0;
%% STEP 1: check data vs feature calculation
% only if the array FEAT is given as an input, this is important
% check number of input arguments > 0

% assign variables of function
if nargin == 3 
    FEAT = varargin{2};     % feature definition is second input argument
    ACT = varargin{1};      % activity data table is first input argument 
elseif nargin == 2 && istable(varargin{1}) == 1
    ACT = varargin{1};      % input is activity data table
elseif nargin == 2 && istable(varargin{1}) == 0
    FEAT = varargin{1};     % input is feature array
end



% if no user-defined sensor selection, base on available data
if exist('FEAT','var') == 0
    FEAT  = [1 2 3 4 5 6 7];        % if no sensor selection provided for feature calculation, base on data available, and thus start with all features
end

% check the availability of data and adjust FEAT accordingly
L = length(DATA.CowID);
% check MILK YIELD data -- less than 80% of measurements have value
if  (isempty(find(string(DATA.Properties.VariableNames) == 'TMY', 1)) == 1 || length(find(isnan(DATA.TMY)==0)) < 0.8*L) && ismember(1,FEAT) == 1
    warning('No TMY variable or not enough data for TMY, TMY features excluded')
    FEAT(FEAT == 1) = [];
end

% check the availability of ACTIVITY data, given the existence of ACT; we
% assume proper preprocessing of the activity data here
if exist('ACT','var') == 0 && ismember(2,FEAT) == 1
    warning('Activity variable not included, activity features excluded')
    FEAT(FEAT == 2) = [];
end

% check WEIGHT data, based on presence of CowWeight var, 60% of data
if  (isempty(find(string(DATA.Properties.VariableNames) == 'CowWeight', 1)) == 1 || length(find(isnan(DATA.CowWeight)==0)) < 0.6*L) && ismember(3,FEAT) == 1
    warning('No WEIGHT variable or not enough data for CowWeight, WEIGHT features excluded')
    FEAT(FEAT == 3) = [];
end

% check ELECTRICAL CONDUCTIVITY data, based on presence of ECLF var, 80%
if  (isempty(find(string(DATA.Properties.VariableNames) == 'ECLF', 1)) == 1 || length(find(isnan(DATA.ECLF)==0)) < 0.8*L) && ismember(4,FEAT) == 1
    warning('No EC variable or not enough data for EC, EC features excluded')
    FEAT(FEAT == 4) = [];
end

% check SCC data, based on presence of SCC var, 30% of data
if  (isempty(find(string(DATA.Properties.VariableNames) == 'SCC', 1)) == 1 || length(find(isnan(DATA.SCC)==0)) < 0.3*L) && ismember(5,FEAT) == 1
    warning('No SCC variable or not enough data for SCC, SCC features excluded')
    FEAT(FEAT == 5) = [];
end

% check FAT/PROT data, based on presence of FAT/PROT var, 30% of data
if  (isempty(find(string(DATA.Properties.VariableNames) == 'Fat', 1)) == 1 || length(find(isnan(DATA.Fat)==0)) < 0.3*L  ...
    || isempty(find(string(DATA.Properties.VariableNames) == 'Protein', 1)) == 1 || length(find(isnan(DATA.Protein)==0)) < 0.3*L) ...
    && ismember(6,FEAT) == 1
    warning('No Fat/Protein variable or not enough data for Fat/Protein, Fat/Protein features excluded')
    FEAT(FEAT == 6) = [];
end

% check TEMPERATURE data, based on presence of MilkT var, 60% of data
if  (isempty(find(string(DATA.Properties.VariableNames) == 'MilkT', 1)) == 1 || length(find(isnan(DATA.MilkT)==0)) < 0.6*L) && ismember(7,FEAT) == 1
    warning('No TEMP variable or not enough data for MilkT, TEMPERATURE features excluded')
    FEAT(FEAT == 7) = [];
end


%% STEP 2: Define CSF and CLSF the features should be added to these tables
% find needed variables
VN1 = find(string(DATA.Properties.VariableNames) == 'CowID');    % col position of CowID
VN2 = find(string(DATA.Properties.VariableNames) == 'Lac');      % col position of Lac
VN3 = find(string(DATA.Properties.VariableNames) == 'Calving');  % col position of Calving

% overview of all available cow lactations, based on CowID and Lac (and
% Calving)
CLSF = unique(DATA(:,[VN1 VN2 VN3]),'rows');   % Summary of unique rows [CowID Lac Calving]
CLSF = sortrows(CLSF,[ 1 3]);               % sort rows of CLSF summary
CLSF.FirstMilk(:,1) = NaN;                    % DIM first milking
CLSF.LastMilk(:,1) = NaN;                     % DIM last milking
CLSF.Start(:,1) = NaT;                        % Startdate for this lactation
CLSF.End(:,1) = NaT;                          % Enddate for this lactation
for i = 1:length(CLSF.Lac)
    ind = find(DATA.Lac == CLSF.Lac(i) & strcmp(DATA.CowID(:,1),CLSF.CowID{i,1}) == 1 &  isnan(DATA.TMY) == 0);
    if isempty(ind) ==0 % if there's data available for this lactation, and TMY not 0 for all
        CLSF.FirstMilk(i,1) = min(DATA.DIM(ind));              % DIM first milking
        CLSF.LastMilk(i,1) = max(DATA.DIM(ind));               % DIM last milking
        CLSF.Start(i,1) = min(DATA.Date(ind));                 % Startdate for this lactation
        CLSF.End(i,1) = max(DATA.Date(ind));                   % Enddate for this lactation
    end
end
clear i ind

% make overview of cows in the data, in analogy with the 'cows' variable to
% containing the rankings. It will be able to merge them in a later stage
% (outside this function).
% overview per cow
CSF = unique(DATA(:,VN1));     % Summary of unique rows [CowID]
CSF.FirstLac(:,1) = NaN;       % First lactation for which data is available
CSF.LastLac(:,1) = NaN;        % Last lactation for which data is available
CSF.Start(:,1) = NaT;          % Start date of measurements for this cow
CSF.End(:,1) = NaT;            % End date of measurments for this cow
CSF.Ndays(:,1) = NaN;          % Number of days data is available for this cow
for i = 1:length(CSF.CowID)
    ind = find(strcmp(CLSF.CowID(:,1),CSF.CowID{i,1}) == 1 & isnan(CLSF.FirstMilk)==0); % Overview based on CLSF, exclude lactations without data
    if isempty(ind) == 0
        CSF.FirstLac(i,1) = min(CLSF.Lac(ind));                   % First lactation
        CSF.LastLac(i,1) = max(CLSF.Lac(ind));                    % Last lactation
        CSF.Start(i,1) = min(CLSF.Start(ind));                    % First measurement of that cow
        CSF.End(i,1) = max(CLSF.End(ind));                        % Last measurment of that cow
        CSF.Ndays(i,1) =  max(datenum(CLSF.End(ind)))-min(datenum(CLSF.Start(ind))); % number of days included
    end
end
clear i ind 

% If no sensor features can be calculated, return function and give warning
if isempty(FEAT) == 1
    warning('No sensor features can be calculated')
    return
end
    
    

%% STEP 3: SENSOR FEATURES MILK YIELD
% if MILK YIELD is included as sensor data, this part of the function will
% calculate sensor features for milk yield characteristics
%   Production performance compared to herd mates in the same lactation
%       in total (actually not sound to be used- in ranking)
%       in first 50 days
%   Lactation curve characteristics (of first 305 days?)
%       slope in first 50 days (WOOD - DIJKSTRA)
%       persistency in linear part after day 65 (WOOD - DIJKSTRA)
%       peak yield  (WOOD - DIJKSTRA)       ==> corrected for average lactation chars?
%       moment of lactation peak (WOOD - DIJKSTRA)
%       peak yield as full maximum of production (RAW DATA)
%       moment of lactation peak (RAW DATA)
%       average RMSE 
% 	Characterization of perturbations
%       N times re-iteration of lactation model is necessary to obtain fit
%       N phases of > 5 days negative residuals
%       N phases of > 10 days negative residuals
%       % negative residuals in first 50 days compared to corrected
%       intercept model
%       Average "loss" of (3) largest periods of residual(s)
%       Number of sign changes of model fit - if high = good fit
%       (comparable with autocorrelation measure) / average autocorrelation
%       of the residuals
%   Characteristics which are related to recalving and dry-off
%       length of the lactation period
%       timing when yield starts decreasing compared to persistent curve
%       length of dry period between two lactations
%   Characteristics of SMOOTHING FUNCTION --> to have a similar rather
%       black box approach for all sensor data

% ADDITIONAL TO INCLUDE LATER/ PREDICTION PERFORMANCE IF USING MM AND BAYES
% PREDICTIONS & comparison of random effects (e.g. using splines MM)- e.g.
% with 5 points - need to see whether there is any correlation in this

if ismember(1, FEAT) == 1           % check whether MY features are calculated
    
    % define 'user defined' functions able to describe the lactation course
    Dijkstra = @(p,t) p(1).*exp((p(2)./p(3))*(1-exp(-p(3).*t))-p(4).*t);
    Wood = @(p,t) p(1).*t.^p(2).*exp(-p(3).*t);
    options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',5000);
    
    % check for DIM variable
    if isempty(find(string(DATA.Properties.VariableNames) == 'DIM', 1)) ==1     % check if DIM variable is available; if no DIM, calculate DIM
        DATA.DIM(:,1) = datenum(DATA.Date(:,1))-datenum(DATA.Calving(:,1)); % prepare array
    end
    
    % set up a structure containing the herd and lactation baselines
    HerdFeat.AverageMilkYield305 = zeros(max(DATA.Lac),306);              % prepare array, row = lac, col = DIM
    HerdFeat.MedianMilkYield305 = zeros(max(DATA.Lac),306);               % prepare array, row = lac, col = DIM
    HerdFeat.WoodLac = zeros(max(DATA.Lac),306);                          % prepare array, row = lac, col = DIM
    HerdFeat.DijkLac = zeros(max(DATA.Lac),306); 
    HerdFeat.WoodPar = zeros(max(DATA.Lac),3);
    HerdFeat.DijkPar = zeros(max(DATA.Lac),4);
    for i = 1:max(DATA.Lac)         % for all lactations in the dataset
        for j = 0:305               % for all days (0 to 305)
            if i < 3
                ind = find(DATA.Lac == i & DATA.DIM == j & DATA.TMY ~= 0);      % separate means for Lac == 1 or Lac == 2
            else
                ind = find(DATA.Lac >= 3 & DATA.DIM == j & DATA.TMY ~= 0);      % from 3 and higher, use average over all Lac >3
            end
            HerdFeat.AverageMilkYield305(i,j+1) = nanmean(DATA.TMY(ind));       % calculate mean for this DIM and lactation
            HerdFeat.MedianMilkYield305(i,j+1) = nanmedian(DATA.TMY(ind));      % calculate mean for this DIM and lactation
        end
        DD = 0:305;
        HerdFeat.WoodPar(i,:) = lsqcurvefit(Wood,[15 0.20 0.004],DD(isnan(HerdFeat.AverageMilkYield305(i,:))==0),HerdFeat.AverageMilkYield305(i,isnan(HerdFeat.AverageMilkYield305(i,:))==0),[0 0 0],[40 1 0.01]);    % Wood parameters
        HerdFeat.WoodLac(i,:) = Wood(HerdFeat.WoodPar(i,:),0:305);              % Wood model values

        HerdFeat.DijkPar(i,:) = lsqcurvefit(Dijkstra,[17 0.08 0.07 0.002],DD(isnan(HerdFeat.AverageMilkYield305(i,:))==0),HerdFeat.AverageMilkYield305(i,isnan(HerdFeat.AverageMilkYield305(i,:))==0),[0 0 0 0],[40 0.5 0.5 0.010]); % Dijkstra pars
        HerdFeat.DijkLac(i,:) = Dijkstra(HerdFeat.DijkPar(i,:),0:305);          % Dijkstra model values
    end
    
    % CAT1 - prepare columns containing the sensor features
    CLSF.TMY305D(:,1) = NaN;        % TMY 305 deviation compared to lactation average, expressed in % of average
    CLSF.TMY50D(:,1) = NaN;         % TMY 50 days deviation compared to lactation average
    CLSF.PeakYieldW(:,1) = NaN;     % Peak MY as calculated max Wood
    CLSF.PeakDIMW(:,1) = NaN;       % DIM Peak MY as calculated max Wood
    CLSF.PeakYieldD(:,1) = NaN;     % Peak MY as calculated max Dijkstra
    CLSF.PeakDIMD(:,1) = NaN;       % DIM Peak MY as calculated max Dijkstra
    CLSF.SlopeInW(:,1) = NaN;       % Slope of MY in first increasing period based on Wood (slope > 2 kg/day)
    CLSF.SlopeInD(:,1) = NaN;       % Slope of MY in first increasing period based on Dijkstra (slope > 2 kg/day)
    CLSF.PersistW(:,1) = NaN;       % Slope of MY in decrease phase
    CLSF.PersistD(:,1) = NaN;       % Slope of MY in decrease phase
    CLSF.MaxYield(:,1) = NaN;       % Peak MY as max raw TMY data
    CLSF.MaxDIM(:,1) = NaN;         % DIM max TMY
    CLSF.RMSEDijk(:,1) = NaN;       % RMSE of Dijkstra model
    CLSF.RMSEWood(:,1) = NaN;       % RMSE of Wood model
    CLSF.pWOOD(:,1:3) = NaN*zeros(length(CLSF.CowID),3);    % Wood pars
    CLSF.pDIJK(:,1:4) = NaN*zeros(length(CLSF.CowID),4);    % Dijkstra pars

    % CAT2A - prepare columns containing the perturbation characterization
    % IN THE FIRST 305 DAYS OF LACTATION
    % USING WOOD
    CLSF.N5dayNegW(:,1) = NaN;      % N phases of > 5 days negative residuals and at least once below 80% of model
    CLSF.N10dayNegW(:,1) = NaN;     % N phases of > 10 days negative residuals and at least once below 80% of model
    CLSF.N1dayNesW(:,1) = NaN;      % percentage of days < 0.85% Wood in first 305days
    CLSF.ARes3per(:,1) = NaN;       % Average "loss" of (3) largest perturbations when production compared to HERD LACTATION average  
    CLSF.ARes3perW(:,1) = NaN;      % Average "loss" of (3) largest perturbations when production compared to WOOD model  
    CLSF.maxResW(:,1) = NaN;        % Absolute max negative residual within period of at least 5 negatives
    CLSF.NiterW(:,1) = NaN;         % N times re-iteration of lactation model is necessary to obtain unchanged fit (deletion of res<80% model; max 20
    CLSF.AutocorrW(:,1) = NaN;      % 'Autocorrelation' expressed as number of sign changes in the residuals
    CLSF.MeanResW(:,1) = NaN;       % Mean of the residuals Wood
    CLSF.AbsMeanResW(:,1) = NaN;    % Mean of the absolute value of the residuals Wood
    CLSF.VarResW(:,1) = NaN;        % Variance of the residuals Wood
    
    % CAT2B - prepare columns containing the perturbation characterization
    % IN THE FIRST 305 DAYS OF LACTATION
    % USING DIJKSTRA
    CLSF.N5dayNegD(:,1) = NaN;      % N phases of > 5 days negative residuals and at least once below 80% of model
    CLSF.N10dayNegD(:,1) = NaN;     % N phases of > 10 days negative residuals and at least once below 80% of model
    CLSF.N1dayNesD(:,1) = NaN;      % percentage of days < 0.85% Dijkstra in first 305days
    CLSF.ARes3perDherd(:,1) = NaN;  % Average "loss" of (3) largest perturbations when production compared to HERD LACTATION average  
    CLSF.ARes3perD(:,1) = NaN;      % Average "loss" of (3) largest perturbations when production compared to WOOD model  
    CLSF.maxResD(:,1) = NaN;        % Absolute max negative residual within period of at least 5 negatives
    CLSF.NiterD(:,1) = NaN;         % N times re-iteration of lactation model is necessary to obtain unchanged fit (deletion of res<80% model; max 20
    CLSF.AutocorrD(:,1) = NaN;      % 'Autocorrelation' expressed as number of sign changes in the residuals
    CLSF.MeanResD(:,1) = NaN;       % Mean of the residuals Dijkstra
    CLSF.AbsMeanResD(:,1) = NaN;    % Mean of the absolute value of the residuals Dijkstra
    CLSF.VarResD(:,1) = NaN;        % Variance of the residuals Dijkstra
   
    
    % CAT3 - characteristics of the perturbation for the 'FULL' lactation 
        % figure, Wood, Dijkstra, but not yet implemented
        
    % CAT4 - Characteristics which are related to recalving and dry-off
    CLSF.LacLength(:,1) = NaN;      % Length of the productive period
    CLSF.StartDecrW(:,1) = NaN;     % Moment that lactation consistently lower than linear part of lactationmodel (iterate!)
    CLSF.StartDecrD(:,1) = NaN;     % Moment that lactation consistently lower than linear part of lactationmodel (iterate!)
    CLSF.LengthDry(:,1) = NaN;      % Based on data availability
    CLSF.EstConc(:,1) = NaN;        % Estimated DIM the cow conceived based on next lac startdate
    CLSF.UNRELI(:,1) = NaN;         % if UNRELI equals 1, this means that the data is suspect in terms of completeness
    
%   Characteristics of SMOOTHING FUNCTION --> to have a similar rather
%       black box approach for all sensor data

    
    % calculate features per lactation
    for i = 1:length(CLSF.CowID)    % for all lactations in CLSF
        ind = find(strcmp(DATA.CowID(:,1),CLSF.CowID{i,1}) == 1 & DATA.Lac(:,1) == CLSF.Lac(i));    % find all data of cow and lac i
%         i;
        if length(ind)<= 4
            CLSF.UNRELI(i,1) = 1;           % data is unreliable bc too little data points
        else
            
            
            % CAT 1 - LACTATION CURVE CHARACTERISTICS
            idx = find(DATA.DIM(ind) < 306);            % all indices for days smaller than 306
            DIMmax = DATA.DIM(ind(idx(end)))+1;         % Days in milk for which data is available
            CLSF.TMY305D(i,1) = -(1-(nansum(DATA.TMY(ind(idx)))/nansum(HerdFeat.AverageMilkYield305(CLSF.Lac(i),1:DIMmax))))*100;

            idx = find(DATA.DIM(ind) < 51);             % all indices for days smaller than 51
            if ~isempty(idx)
                DIMmax = DATA.DIM(ind(idx(end)))+1;         % Days in milk for which data is available
                CLSF.TMY50D(i,1) = -(1-(nansum(DATA.TMY(ind(idx)))/nansum(HerdFeat.AverageMilkYield305(CLSF.Lac(i),1:DIMmax))))*100;
            end
            
            idx = find(DATA.DIM(ind) < 306 & isnan(DATA.TMY(ind))==0);            % all indices for days smaller than 306
            TMY = DATA.TMY(ind(idx));                    % TMY used for modelling
            DIM = DATA.DIM(ind(idx));                    % DIM used for modelling
            pDIJK = lsqcurvefit(Dijkstra,[15 0.04 0.04 0.001],DIM,TMY,[0 0 0 0],[100 10 10 10],options);
            pWOOD = lsqcurvefit(Wood,[15 0.20 0.004],DIM,TMY,[0 0 0],[40 1 0.01],options);

            fWOOD = Wood(pWOOD,DIM);                     % function values of Wood
            fDIJK = Dijkstra(pDIJK,DIM);                 % function values of Dijkstra

            [CLSF.PeakYieldW(i,1), idx] = max(fWOOD);    % Peak MY as calculated max WOOD
            CLSF.PeakDIMW(i,1) = DIM(idx);               % DIM Peak MY as calculated max Dijkstra

            [CLSF.PeakYieldD(i,1), idx] = max(fDIJK);    % Peak MY as calculated max WOOD
            CLSF.PeakDIMD(i,1) = DIM(idx);               % DIM Peak MY as calculated max Dijkstra

% % % % %             F1 = figure('Units','centimeters','OuterPosition',[2 2 26 14]); % subplot(2,2,1); 
% % % % %                         plot(DIM, TMY,'o-','MarkerSize',6, 'LineWidth',2,'Color',[0.73 0.83 0.96]); hold on
% % % % %     %                     plot(DIM,Dijkstra([25.6 0.0975 0.0006 0.09], DIM),'--','LineWidth',1.2)
% % % % %                         plot(DIM,Dijkstra(pDIJK,DIM),'LineWidth',3, 'Color',[0.08 0.17 0.55]); 
% % % % %     %                     plot(DIM,Wood([15 0.20 0.004],DIM),'--','LineWidth',1.2);
% % % % %                         plot(DIM,Wood(pWOOD,DIM),'LineWidth',2.5, 'Color',[48/255 184/255 48/255]);
% % % % %     %                     legend({'TMY','Dijkstra_i_n_i','Dijkstra_f_i_t','Wood_i_n_i','Wood_f_i_t'},'Location','best'); 
% % % % %                         xlabel('DIM (days)'); ylabel('TMY (kg)'); xlim([0 306])
% % % % %                         title(['Lactation features for cow ' CLSF.CowID{i,1} ' in lac ' num2str(CLSF.Lac(i))])
% % % % %                         plot(CLSF.PeakDIMW(i,1),CLSF.PeakYieldW(i,1),'o','LineWidth',2.5,'MarkerSize',8, 'Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0])
% % % % %                         plot(CLSF.PeakDIMD(i,1),CLSF.PeakYieldD(i,1),'o','LineWidth',2.5,'MarkerSize',8, 'Color','b','MarkerFaceColor','b')    

            dWOOD = diff(fWOOD);                         % find slope of Wood model
            idx = find(dWOOD > 0.5,1,'last');            % last measurement with slope > 0.5 kg/day
            if isempty(idx)==0
                CLSF.SlopeInW(i,1) = (fWOOD(idx)-fWOOD(1))/(DIM(idx)-DIM(1));    % Slope of MY in increasing part
% % % % %                         plot([DIM(1) DIM(idx)],[fWOOD(1) fWOOD(idx)],'-.','LineWidth',2.5, 'Color',[0.47 0.67 0.19])
            end            
            
            dDIJK = diff(fDIJK);  
            idx = find(dDIJK > 0.5,1,'last');              % last measurement with slope > 0.5 kg/day
            if isempty(idx) ==0
                CLSF.SlopeInD(i,1) = (fDIJK(idx)-fDIJK(1))/(DIM(idx)-DIM(1));    % Slope of MY in increasing part
% % % % %                         plot([DIM(1) DIM(idx)],[fDIJK(1) fDIJK(idx)],'b-.','LineWidth',2.5)
            end

            idx = find(dWOOD < 0.05,1,'first');          % first measurement with rather flat slope
            if isempty(idx) ==0
                CLSF.PersistW(i,1) = (fWOOD(end)-fWOOD(idx))/(DIM(end)-DIM(idx));    % Slope of MY in day 65 - 305
% % % % %                         plot([DIM(idx) DIM(end)],[fWOOD(idx) fWOOD(end)],'-.','LineWidth',2.5, 'Color',[0.47 0.67 0.19])
            end

                          % find slope of Wood model
            idx = find(dDIJK < 0.05,1,'first');         % first measurement with rather flat slope
            if isempty(idx) ==0
                CLSF.PersistD(i,1) = (fDIJK(end)-fDIJK(idx))/(DIM(end)-DIM(idx));    % Slope of MY in day 65 - 305
% % % % %                         plot([DIM(idx) DIM(end)],[fDIJK(idx) fDIJK(end)],'b-.','LineWidth',2.5)
            end

            [CLSF.MaxYield(i,1),idx] = max(TMY);        % Peak MY as max raw TMY data
            CLSF.MaxDIM(i,1) = DIM(idx);                % DIM max TMY
% % % % %                         plot(CLSF.MaxDIM(i,1),CLSF.MaxYield(i,1),'o','LineWidth',2.5,'MarkerSize',8,'MarkerFaceColor',[0.64 0.08 0.18],'Color',[0.64 0.08 0.18])    

            CLSF.RMSEDijk(i,1) = sqrt(mean((TMY-fDIJK).^2));  % RMSE of Dijkstra model
            CLSF.RMSEWood(i,1) = sqrt(mean((TMY-fWOOD).^2));  % RMSE of Wood Model

            CLSF.pWOOD(i,1:3) = pWOOD;                  % Wood pars
            CLSF.pDIJK(i,1:4) = pDIJK;                  % Dijkstra pars

% % % % %             legend({'TMY','Dijkstra_f_i_t','Wood_f_i_t','Peak_W','Peak_D','Slope_5_0_W','Slope_5_0_D','Slope_l_i_n_W','Slope_l_i_n_D','Max Yield'},'Location','best','FontSize',8);

            % CATEGORY 2A - CHARACTERIZATION OF PERTURBATIONS IN THE FIRST 305 DAYS - WOOD BASED
% % % % %             F2 = figure('Units','centimeters','OuterPosition',[1 1 26 18]);  
% % % % %                         subplot(2,1,1); hold on
% % % % %                         plot(DIM, TMY,'o-','MarkerSize',6, 'LineWidth',2,'Color',[0.73 0.83 0.96]); 
% % % % %                         plot(DIM,Wood(pWOOD,DIM),'LineWidth',2.5, 'Color',[48/255 184/255 48/255]);
% % % % %     %                     legend({'TMY','Dijkstra_i_n_i','Dijkstra_f_i_t','Wood_i_n_i','Wood_f_i_t'},'Location','best'); 
% % % % %                         xlabel('DIM (days)'); ylabel('TMY (kg)'); xlim([0 max(306,max(DIM)+1)])
% % % % %                         title(['Yield perturbation features Wood for cow ' CLSF.CowID{i,1} ' in lac ' num2str(CLSF.Lac(i))])

            RWOOD = TMY - fWOOD;            % Residuals of the Wood model
% % % % %                         subplot(2,1,2); hold on; xlabel('DIM (days)'); ylabel('TMY residual (kg)'); xlim([0 305])
% % % % %                         title(['Residual yield perturbation features Wood for cow ' CLSF.CowID{i,1} ' in lac ' num2str(CLSF.Lac(i))])
% % % % % 
% % % % %                         plot(DIM, zeros(length(DIM),1),'k--', 'LineWidth',1.2); 
% % % % %                         plot(DIM, RWOOD,'o-','MarkerSize',6, 'LineWidth',2,'Color',[0.39 0.47 0.64]); 

            idx = find(TMY < 0.85*fWOOD); 
% % % % %                         subplot(2,1,1); plot(DIM(idx),TMY(idx),'ro','LineWidth',2,'MarkerSize',6)
% % % % %                         subplot(2,1,2); plot(DIM(idx),RWOOD(idx),'ro','LineWidth',2,'MarkerSize',6)

            NegRes = (RWOOD<0);             % 1 for all negative Wood residuals

            [~, idx2] = F3_pattern(NegRes,[1 1 1 1 1]);   % indices for which 5 consecutive measurements the residual is negative
            Negs = zeros(length(idx2),6);
            for j = 1:length(idx2)
                idx3 = find(DIM>DIM(idx2(j)) & RWOOD > 0,1,'first');
                if isempty(idx3); idx3 = length(RWOOD);end
                if idx3 ~= length(RWOOD)
                    Negs(j,1:2) = [DIM(idx2(j)) DIM(idx3)-1 ];
                else
                    Negs(j,1:2) = [DIM(idx2(j)) DIM(idx3)];
                end

                Negs(j,3) = max(ismember(idx2(j):idx3,idx));
                Negs(j,4) = DIM(idx3)-1-DIM(idx2(j));
% % % % %                         if Negs(j,3) == 1 && Negs(j,4) <10
% % % % %                             subplot(2,1,1); plot([Negs(j,1) Negs(j,1)],[0 60],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[0 60],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                             subplot(2,1,2); plot([Negs(j,1) Negs(j,1)],[-20 20],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[-20 20],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                         elseif Negs(j,3) == 1 && Negs(j,4) >= 10
% % % % %                             subplot(2,1,1); plot([Negs(j,1) Negs(j,1)],[0 60],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[0 60],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                             subplot(2,1,2); plot([Negs(j,1) Negs(j,1)],[-20 20],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[-20 20],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                         end
                 if idx2(j)-idx3 ~= DIM(idx3)-DIM(idx2(j))
                     Negs(j,5) = NaN;
                 else
                     Negs(j,5) = sum(TMY(idx2(j):idx3)-HerdFeat.AverageMilkYield305(CLSF.Lac(i),DIM(idx2(j)):DIM(idx3))');
                 end
                 Negs(j,6) = sum(TMY(idx2(j):idx3)-fWOOD(idx2(j):idx3));    % 
            end
            CLSF.N5dayNegW(i,1) = nansum(Negs(:,3));             % N phases of > 5 days negative residuals and at least once below 85% of model value
            CLSF.N10dayNegW(i,1) = sum(Negs(Negs(:,4)>=10,3));   % N phases of > 10 days negative residuals and at least once below 85% of model
            CLSF.N1dayNesW(i,1) = (length(idx)/length(DIM))*100; % percentage of residuals below 85% of wood compared to N of measurments
            
% % % % %                         subplot(2,1,1); plot(0:DIM(end),HerdFeat.AverageMilkYield305(CLSF.Lac(i),1:DIM(end)+1),':','LineWidth',2,'Color',[0.8 0.8 0.8]);

            Negs = sortrows(Negs, 5,'ascend');                  % sort by 'size' of perturbation in terms of total loss compared to mean
            CLSF.ARes3per(i,1) = sum(Negs(1:min(3,j),5))/sum(Negs(1:min(3,j),4));     % Average "loss" of (3) largest perturbations, when production compared to HERD    

            Negs = sortrows(Negs, 6,'ascend');                  % sort by 'size' of perturbation in terms of total loss compared WOOD
            CLSF.ARes3perW(i,1) = sum(Negs(1:min(3,j),6))/sum(Negs(1:min(3,j),4));    % Average "loss" of (3) largest perturbations when production compared to WOOD model  

            CLSF.maxResW(i,1) = min(RWOOD(4:end));              % minimal residual

            CLSF.MeanResW(i,1) = nanmean(RWOOD);                % Mean of the residuals Wood
            CLSF.AbsMeanResW(i,1) = nanmean(abs(RWOOD));        % Mean of the absolute value of the residuals Wood
            CLSF.VarResW(i,1) = nanvar(RWOOD);                  % Variance of the residuals Wood

            
            DIM2 = DIM;     % prepare DIM
            TMY2 = TMY;     % prepare TMY
            fWOOD2 = fWOOD; % prepare fWOOD
            RMSE = sqrt(mean((TMY2-fWOOD2).^2));        % RMSE
            RMSE2 = 0;
            n=0;
            while abs(RMSE-RMSE2) > 0.1  && n <= 20     % while the difference between the model and the previous model is larger than 0.1 (??)
                n=n+1;                                  % and also, at most 20 iterations
                DIM2 = DIM2(TMY2>0.85*fWOOD2);          % DIM of all milkings above 0.85 of the prediction
                TMY2 = TMY2(TMY2>0.85*fWOOD2);          % TMY of all milkings above 0.85 of the prediction

                if length(DIM2)>3
                    p = lsqcurvefit(Wood,pWOOD,DIM2,TMY2,[0 0 0],[40 1 0.01],options);  % refit model
                    fWOOD2 = Wood(p,DIM2);                  % calculate model values

                    RMSE2 = RMSE;                           % Set the RMSE to the previous ones
                    RMSE = sqrt(mean((TMY2-fWOOD2).^2));    % RMSE
                else
                    n = 100;
                end
            end

% % % % %                         subplot(2,1,1); plot(DIM2,TMY2,'ko--','MarkerSize',2,'LineWidth',1.2)
% % % % %                         plot(DIM,Wood(p,DIM),'-','LineWidth',1.4)

            CLSF.NiterW(i,1) = n-1;         % N times re-iteration of lactation model is necessary to obtain unchanged fit (deletion of res<80% model; max 20

            CLSF.AutocorrW(i,1) = (1-length(find(diff(sign(RWOOD))>1))/length(ind))*100;       % 'Autocorrelation' expressed as number of sign changes in the residuals

            % CATEGORY 2B - CHARACTERIZATION OF PERTURBATIONS IN THE FIRST 305 DAYS - DIJKSTRA BASED
% % % % %             F3 = figure('Units','centimeters','OuterPosition',[2 2 26 14]);  subplot(2,1,1); hold on
% % % % %                         plot(DIM, TMY,'o-','MarkerSize',6, 'LineWidth',2,'Color',[0.73 0.83 0.96]); 
% % % % %                         plot(DIM,Dijkstra(pDIJK,DIM),'LineWidth',3, 'Color',[0.08 0.17 0.55]); 
% % % % %                         xlabel('DIM (days)'); ylabel('TMY (kg)'); xlim([0 max(306,max(DIM)+1)])
% % % % %                         title(['Yield perturbation features Dijkstra for cow ' CLSF.CowID{i,1} ' in lac ' num2str(CLSF.Lac(i))])

            RDIJK = TMY - fDIJK;            % Residuals of the Dijkstra model
% % % % %                         subplot(2,1,2); hold on; xlabel('DIM (days)'); ylabel('TMY residual (kg)'); xlim([0 305])
% % % % %                         title(['Residual yield perturbation features Dijkstra for cow ' CLSF.CowID{i,1} ' in lac ' num2str(CLSF.Lac(i))])
% % % % % 
% % % % %                         plot(DIM, zeros(length(DIM),1),'k--', 'LineWidth',1.2); 
% % % % %                         plot(DIM, RDIJK,'o-','MarkerSize',6, 'LineWidth',2,'Color',[0.39 0.47 0.64]); 

            idx = find(TMY < 0.85*fDIJK); 
% % % % %                         subplot(2,1,1); plot(DIM(idx),TMY(idx),'ro','LineWidth',2,'MarkerSize',6)
% % % % %                         subplot(2,1,2); plot(DIM(idx),RDIJK(idx),'ro','LineWidth',2,'MarkerSize',6)

            NegRes = (RDIJK<0);             % 1 for all negative Wood residuals

            [~, idx2] = F3_pattern(NegRes,[1 1 1 1 1]);   % indices for which 5 consecutive measurements the residual is negative
            Negs = zeros(length(idx2),6);                 % prepare array for negative residuals
            for j = 1:length(idx2)
                idx3 = find(DIM>DIM(idx2(j)) & RDIJK > 0,1,'first'); % find the length of the negative perturbation
                if isempty(idx3); idx3 = length(RDIJK);end           % if the perturbation lasts to end, end = last index
                if idx3 ~= length(RDIJK)
                    Negs(j,1:2) = [DIM(idx2(j)) DIM(idx3)-1 ];       % fill in last negative residual
                else
                    Negs(j,1:2) = [DIM(idx2(j)) DIM(idx3)];
                end

                Negs(j,3) = max(ismember(idx2(j):idx3,idx));         % find whether the 'perturbation' contains a neg < 0.85*model
                Negs(j,4) = DIM(idx3)-1-DIM(idx2(j));                % fill in length of the perturbation
% % % % %                         if Negs(j,3) == 1 && Negs(j,4) <10
% % % % %                             subplot(2,1,1); plot([Negs(j,1) Negs(j,1)],[0 60],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[0 60],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                             subplot(2,1,2); plot([Negs(j,1) Negs(j,1)],[-20 20],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[-20 20],'--','LineWidth',2,'Color',[0.51 0.38 0.48])
% % % % %                         elseif Negs(j,3) == 1 && Negs(j,4) >= 10
% % % % %                             subplot(2,1,1); plot([Negs(j,1) Negs(j,1)],[0 60],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[0 60],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                             subplot(2,1,2); plot([Negs(j,1) Negs(j,1)],[-20 20],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                                             plot([Negs(j,2) Negs(j,2)],[-20 20],'--','LineWidth',2,'Color',[0.58 0.39 0.39])
% % % % %                         end
                 if idx2(j)-idx3 ~= DIM(idx3)-DIM(idx2(j))
                     Negs(j,5) = NaN;
                 else
                     Negs(j,5) = sum(TMY(idx2(j):idx3)-HerdFeat.AverageMilkYield305(CLSF.Lac(i),DIM(idx2(j)):DIM(idx3))');
                 end
                 Negs(j,6) = sum(TMY(idx2(j):idx3)-fDIJK(idx2(j):idx3));    % fill in average loss compared to DIJKSTRA model
            end
            CLSF.N5dayNegD(i,1) = nansum(Negs(:,3));                % N phases of > 5 days negative residuals and at least once below 90% of model value
            CLSF.N10dayNegD(i,1) = nansum(Negs(Negs(:,4)>=10,3));   % N phases of > 10 days negative residuals and at least once below 80% of model
            CLSF.N1dayNesD(i,1) = (length(idx)/length(DIM))*100; % percentage of residuals below 85% of wood compared to N of measurments

% % % % %                         subplot(2,1,1); plot(0:DIM(end),HerdFeat.AverageMilkYield305(CLSF.Lac(i),1:DIM(end)+1),':','LineWidth',2,'Color',[0.8 0.8 0.8]);

            Negs = sortrows(Negs, 5,'ascend');                  % sort by 'size' of perturbation in terms of total loss compared to mean
            CLSF.ARes3perD(i,1) = nansum(Negs(1:min(3,j),5))/sum(Negs(1:min(3,j),4));     % Average "loss" of (3) largest perturbations, when production compared to HERD    

            Negs = sortrows(Negs, 6,'ascend');                  % sort by 'size' of perturbation in terms of total loss compared WOOD
            CLSF.ARes3perD(i,1) = nansum(Negs(1:min(3,j),6))/sum(Negs(1:min(3,j),4));    % Average "loss" of (3) largest perturbations when production compared to WOOD model  

            CLSF.maxResD(i,1) = min(RDIJK);     % minimal residual
            CLSF.MeanResD(i,1) = nanmean(RDIJK);                % Mean of the residuals Wood
            CLSF.AbsMeanResD(i,1) = nanmean(abs(RDIJK));        % Mean of the absolute value of the residuals Wood
            CLSF.VarResD(i,1) = nanvar(RDIJK);                  % Variance of the residuals Wood

            DIM2 = DIM;     % prepare DIM
            TMY2 = TMY;     % prepare TMY
            fDIJK2 = fDIJK; % prepare fWOOD
            RMSE = sqrt(mean((TMY2-fDIJK2).^2));        % RMSE
            RMSE2 = 0;
            n=0;
            while abs(RMSE-RMSE2) > 0.1  && n <= 20     % while the difference between the model and the previous model is larger than 0.1 (??)
                n=n+1;                                  % and also, at most 20 iterations
                DIM2 = DIM2(TMY2>0.85*fDIJK2);          % DIM of all milkings above 0.85 of the prediction
                TMY2 = TMY2(TMY2>0.85*fDIJK2);          % TMY of all milkings above 0.85 of the prediction
                if length(TMY2)>5
                    p = lsqcurvefit(Dijkstra,pDIJK,DIM2,TMY2,[0 0 0 0],[100 10 10 10],options);  % refit model
                    fDIJK2 = Dijkstra(p,DIM2);                  % calculate model values

                    RMSE2 = RMSE;                           % Set the RMSE to the previous ones
                    RMSE = sqrt(mean((TMY2-fDIJK2).^2));    % RMSE
                else
                    n = 100;
                end
            end
i
% % % % %                         subplot(2,1,1); plot(DIM2,TMY2,'ko--','MarkerSize',2,'LineWidth',1.2)
% % % % %                         plot(DIM,Dijkstra(p,DIM),'-','LineWidth',1.4)

            CLSF.NiterD(i,1) = n-1;         % N times re-iteration of lactation model is necessary to obtain unchanged fit (deletion of res<80% model; max 20

            CLSF.AutocorrD(i,1) = (1-length(find(diff(sign(RDIJK))>1))/length(ind))*100;       % 'Autocorrelation' expressed as number of sign changes in the residuals


            % CATEGORY 3 - CHARACTERIZATION OF PERTURBATIONS IN THE WHOLE LACTATION PERIOD
            ind = find(strcmp(DATA.CowID(:,1),CLSF.CowID{i,1}) == 1 & DATA.Lac(:,1) == CLSF.Lac(i) & isnan(DATA.TMY)==0);    % find all data of cow and lac i
            DIM = DATA.DIM(ind);                        % all DIM data
            TMY = DATA.TMY(ind);                        % all TMY data
            if DATA.DIM(ind(end)) > 305

                pWOOD = lsqcurvefit(Wood,pWOOD,DIM,TMY,[0 0 0],[40 1 0.01],options);            % fit wood
                fWOOD = Wood(pWOOD,DIM); 
                pDIJK = lsqcurvefit(Dijkstra,pDIJK,DIM,TMY,[0 0 0 0],[100 10 10 10],options);   % fit dijkstra
                fDIJK = Dijkstra(pDIJK,DIM);
            end           
% % % % %             F4 =  figure('Units','centimeters','OuterPosition',[2 2 26 14]); % subplot(2,2,3); 
% % % % %                         plot(DIM, TMY,'o-','MarkerSize',6, 'LineWidth',2,'Color',[0.73 0.83 0.96]); hold on
% % % % %                         plot(DIM,fDIJK,'LineWidth',3, 'Color',[0.08 0.17 0.55]); 
% % % % %                         plot(DIM,fWOOD,'LineWidth',2.5, 'Color',[48/255 184/255 48/255]);
% % % % %                         xlabel('DIM (days)'); ylabel('TMY (kg)'); xlim([0 max(306,max(DIM)+1)])
% % % % %                         title(['General lactation characteristics for cow ' CLSF.CowID{i,1} ' in lac ' num2str(CLSF.Lac(i))])

            % CATEGORY 4 - GENERAL LACTATION CHARACTERISTICS - RECALVING and DRY-OFF
            CLSF.LacLength(i,1) = max(DIM);      % Length of the productive period

            % Wood - based find consistently low
            DIM2 = DATA.DIM(ind);     % prepare DIM - all data
            TMY2 = DATA.TMY(ind);     % prepare TMY - all data
            fWOOD2 = Wood(pWOOD,DIM2); % prepare fWOOD
            RMSE = sqrt(mean((TMY2-fWOOD2).^2));        % RMSE
            RMSE2 = 0;
            n=0;
            while abs(RMSE-RMSE2) > 0.1  && n <= 20     % while the difference between the model and the previous model is larger than 0.1 (??)
                n=n+1;                                  % and also, at most 20 iterations
                DIM2 = DIM2(TMY2>0.85*fWOOD2);          % DIM of all milkings above 0.85 of the prediction
                TMY2 = TMY2(TMY2>0.85*fWOOD2);          % TMY of all milkings above 0.85 of the prediction
                
                if length(DIM2)>4
                    p = lsqcurvefit(Wood,pWOOD,DIM2,TMY2,[0 0 0],[40 1 0.01],options);  % refit model
                    fWOOD2 = Wood(p,DIM2);                  % calculate model values

                    RMSE2 = RMSE;                           % Set the RMSE to the previous ones
                    RMSE = sqrt(mean((TMY2-fWOOD2).^2));    % RMSE
                else
                    n = 100;
                end
            end
            CLSF.StartDecrW(i,1) = DIM2(end);% Start moment that lactation consistently lower than linear part of lactationmodel
% % % % %                         plot(DIM2,TMY2,'ko--','MarkerSize',2,'LineWidth',1.2)
% % % % %                         plot(DIM,Wood(p,DIM),'-','LineWidth',1.4)
% % % % %                         plot([DIM2(end) DIM2(end)],[0 60],'r--','LineWidth',2)

            % Dijkstra -based find consistently low
            DIM2 = DATA.DIM(ind);     % prepare DIM - all data
            TMY2 = DATA.TMY(ind);     % prepare TMY - all data
            fDIJK2 = fDIJK; % prepare fWOOD
            RMSE = sqrt(mean((TMY2-fDIJK2).^2));        % RMSE
            RMSE2 = 0;
            n=0;
            while abs(RMSE-RMSE2) > 0.1  && n <= 20     % while the difference between the model and the previous model is larger than 0.1 (??)
                n=n+1;                                  % and also, at most 20 iterations
                DIM2 = DIM2(TMY2>0.85*fDIJK2);          % DIM of all milkings above 0.85 of the prediction
                TMY2 = TMY2(TMY2>0.85*fDIJK2);          % TMY of all milkings above 0.85 of the prediction
                if length(DIM2)>4
                    p = lsqcurvefit(Dijkstra,pDIJK,DIM2,TMY2,[0 0 0 0],[100 10 10 10],options);  % refit model
                    fDIJK2 = Dijkstra(p,DIM2);                  % calculate model values

                    RMSE2 = RMSE;                           % Set the RMSE to the previous ones
                    RMSE = sqrt(mean((TMY2-fDIJK2).^2));    % RMSE
                else
                    n = 100;
                end
            end
            CLSF.StartDecrD(i,1) = DIM2(end);% Start moment that lactation consistently lower than linear part of lactationmodel
% % % % %                         plot(DIM2,TMY2,'o--','MarkerSize',2,'LineWidth',1.2,'Color', [0.5 0.5 0.5])
% % % % %                         plot(DIM,Dijkstra(p,DIM),'-','LineWidth',1.4)
% % % % %                         plot([DIM2(end) DIM2(end)],[0 60],'r:','LineWidth',2)


            verification = sortrows(CLSF,{'CowID','Lac'});  % make sure that the cows are in the right order in CSLF
            idx = find(strcmp(verification.CowID(:,1),CLSF.CowID{i,1}) == 1);
            if CLSF.Lac(i) ~= max(verification.Lac(idx))
                k = find(verification.Lac(idx) == CLSF.Lac(i));
                CLSF.LengthDry(i,1) = datenum(verification.Start(idx(k+1)))-datenum(CLSF.End(i));      % Based on data availability        
                CLSF.EstConc(i,1) = (datenum(verification.Start(idx(k+1)))-283)-datenum(CLSF.Calving(i)); % Estimated DIM the cow conceived based on next lac startdate
            end
            
            if DIM(end) < 200 && TMY(end) > 8
                CLSF.UNRELI(i,1) = 1;           % data is unreliable bc too little data points
            else 
                CLSF.UNRELI(i,1) = 0;           % data is unreliable bc too little data points
            end

            
% % % % %             saveas(F1,['path\'  'LacFeat_1_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.jpg'],'jpg')
% % % % %             saveas(F1,['path\'  'MLacFeat_1_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.fig'],'fig')
% % % % %             saveas(F2,['path\'  'LacFeat_2_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.jpg'],'jpg')
% % % % %             saveas(F2,['path\'  'MLacFeat_2_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.fig'],'fig')
% % % % %             saveas(F3,['path\'  'LacFeat_3_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.jpg'],'jpg')
% % % % %             saveas(F3,['path\'  'MLacFeatm_3_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.fig'],'fig')
% % % % %             saveas(F4,['path\'  'LacFeat_4_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.jpg'],'jpg')
% % % % %             saveas(F4,['path\'  'MLacFeatm_4_' CLSF.CowID(i) 'Lac_' num2str(CLSF.Lac(i)) '.fig'],'fig')
        end
    end
end


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

if ismember(2, FEAT) == 1   % check whether MY features need to be calculated 
    % predefine the different cols for filling in feature characteristics
    % CAT1 - general activity data charactistics
    CLSF.ACStart(:,1) = NaN;            % Start DIM of the activity data
    CLSF.ACEnd(:,1) = NaN;              % End DIM of the activity data
    CLSF.ACDailyMean(:,1) = NaN;        % Mean total activity for that cow
    CLSF.ACNumberZero(:,1) = NaN;       % Number of days the activity sum is zero
    CLSF.ACPercZero(:,1) = NaN;         % Percentage of zero activity 
    CLSF.ACSkewness(:,1) = NaN;         % Skewness of the activity data
    CLSF.ACOverallVariance(:,1) = NaN;  % Average within-day variance in activity for that cow
    CLSF.ACOverallActHerd(:,1) = NaN;   % The mean activity of that cow compared to her herd mates
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
    
    for i = 1:length(CLSF.CowID)    % for all the cows, find and summarize the activity data
        ind = find(strcmp(ACT.CowID(:,1),CLSF.CowID{i,1}) == 1 & datenum(ACT.Date(:,1))>= datenum(CLSF.Start(i)) & datenum(ACT.Date(:,1))<= datenum(CLSF.End(i)));    % find all data of cow and lac i

        % first step = set up new dataset which 'summarizes' the different
        % activity measures from the raw data
        if isempty(ind) == 0    % check whether we have activity data available for this cow lactation
            CowAct = ACT(ind,:);            % all activity data of this cow
            CowAct.DIM = datenum(CowAct.Date)-datenum(CLSF.Calving(i))+0.0001;     % add the days in milk to the data
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
end  

ALLCAS = array2table(ALLCAS,'VariableNames',{'Actsum','Actmean','Actstd','Actvar','Actmin','Actmax','DIM'});
ALLCAS(ALLCAS.DIM == 0 & ALLCAS.Actsum == 0 & ALLCAS.Actstd == 0 & ALLCAS.Actmax == 0,:) = [];      % delete empty rows
