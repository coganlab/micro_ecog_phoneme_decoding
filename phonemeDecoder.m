addpath(genpath('C:\Users\sd355\Box\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\Data\Micro\Processed Data';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window


baseTW = [-0.5 0]; % preonset time window
activeTW = [-1 1]; % postonset time window
gammaF = [70 150]; % frequency in Hz
fsDown = 200; % sampling frequency in Hz

%% data Loading
subjectId = 'S23';

% Loading Experiment - Experiment contains information necessary to task
% information
Experiment = loadExperiment(subjectId);
% Loading channel Map (2D channel indices)
chanMap = Experiment.recording.channel_map;
chanMap = (chanMap');
selectedChannels = sort(chanMap(~isnan(chanMap)))';
% sampling rate
fsD = Experiment.processing.ieeg.sample_rate;
% Loading trial-ids for the subject and corresponding day
Trials = dbTrials(subjectId,Experiment.recording.recording_day,'Speech_OvertMimeMove');

trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
% Loading data corresponding to phoneme 1 onset & auditory onset
[ieegResponse,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,'phon1Onset','ieeg',tw.*1000);
[ieegAuditory,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory','ieeg',tw.*1000);

% Changing the dimensions from trials x channels x timepoints to channels x
% trials x timepoints
ieegAuditory = permute(ieegAuditory,[2,1,3]);
ieegResponse = permute(ieegResponse,[2,1,3]);
 
% removing non-response data
respId = find(~isnan(trigOnset));
ieegResponse = ieegResponse(:,respId,:);
ieegAuditory = ieegAuditory(:,respId,:); 
Trials = Trials(respId);
trigRespOfs = [Trials.ResponseOffset]./30e3;  
trigOns = trigOnset(respId)./fsD;
trigAud = trigAuditory(respId)./fsD;
respTime  = trigRespOfs-trigOns;
% [respTimeSort,sortId] = sort(respTime);
%% Ieeg Class definition

ieegBase = ieegStructClass(ieegAuditory, fsD, tw, [1 fsD/2], 'Auditory');
ieegResponse = ieegStructClass(ieegResponse, fsD, tw, [1 fsD/2], 'Response');

%% Phoneme trial parser
phonemeTrial = phonemeSequenceTrialParser(Trials);
%% Common average referencing
load([subjectId '_impedance.mat'])
higImpId = find(log10(impedance)>6);
ieegBaseCar = ieegBase.extractCar(higImpId);
clear ieegBase
ieegResponseCar = ieegResponse.extractCar(higImpId);
clear ieegResponse
%% Extract High Gamma
normType =1 ;
ieegHGBase = ieegBaseCar.extractHiGamma(fsDown,baseTW);
normFactor = extractHGnormFactor(ieegHGBase);
ieegHGBaseNorm = normHiGamma(ieegHGBase,normFactor,2);
ieegHGRespNorm = ieegResponseCar.extractHiGamma(fsDown,activeTW,normFactor, normType);
%% Tsne score
load([subjectId '_sigChannel.mat'])
numFold = 20;
varExplained = 80;
chanMapSig = chanMap;
chanMapSig(~ismember(chanMap(:),sigChannel))= nan;
ieegHGRespNormSig = ieegHGRespNorm;
ieegHGRespNormSig.data = ieegHGRespNorm.data(sigChannel,:,:);
decoderObject = decoderClass(numFold,varExplained);
decompStructPhoneme = tsneDecompose(decoderObject,ieegHGRespNormSig,phonemeTrial.phonemeUnit(:,1)',chanMapSig, [-0.25 0.25],[256 128 64 32],50);
decompStructPhonemeChance = tsneDecompose(decoderObject,ieegHGRespNormSig,shuffle(phonemeTrial.phonemeUnit(:,1)'),chanMapSig, [-0.25 0.25],[1024],50);
decompStructClass = tsneDecompose(decoderObject,ieegHGRespNormSig,phonemeTrial.phonemeClass(:,1)',chanMapSig, [-0.25 0.25],[256 128 64 32],50);
decompStructClassChance = tsneDecompose(decoderObject,ieegHGRespNormSig,shuffle(phonemeTrial.phonemeClass(:,1)'),chanMapSig, [-0.25 0.25],[256 128 64 32],50);

%%
w95 = 0.7193;
figure;
subplot(2,1,1);
h = boxplot(decompStructClass.silScore', 'symbol','','Whisker',w95);
hold on;
yline(mean(decompStructClassChance.silScore(:)), '--','chance','LineWidth',1);
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
ylabel('Silhouette score');
title('Articulator grouping');
ylim([0 max(decompStructClass.silScore(:))]);

axis square
set( gca                       , ...
    'FontName'   , 'Arial' );
legend off
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );

subplot(2,1,2);
h = boxplot(decompStructPhoneme.silScore', 'symbol','','Whisker',w95);
hold on;
yline(mean(decompStructPhonemeChance.silScore), '--','chance','LineWidth',1);
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
xlabel('% Sampling');
title('Phoneme grouping');
ylim([0 max(decompStructPhoneme.silScore(:))]);

axis square

set( gca                       , ...
    'FontName'   , 'Arial' );

legend off
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );

figure;
subplot(2,1,1);
h = boxplot(decompStructClass.silRatio', 'symbol','','Whisker',w95);
hold on;
yline(mean(decompStructClassChance.silRatio(:)), '--','chance','LineWidth',1);
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
ylabel('Silhouette ratio');
title('Articulator grouping');
ylim([0 max(decompStructClass.silRatio(:))]);
set( gca                       , ...
    'FontName'   , 'Arial' );

legend off
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );
subplot(2,1,2);
h = boxplot(decompStructPhoneme.silRatio', 'symbol','','Whisker',w95);
hold on;
yline(mean(decompStructPhonemeChance.silRatio), '--','chance','LineWidth',1);
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
xlabel('% Sampling');
title('Phoneme grouping');
ylim([0 max(decompStructPhoneme.silRatio(:))]);

axis square

set( gca                       , ...
    'FontName'   , 'Arial' );

legend off
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'LineWidth'   , 1         );
%% Phoneme decoding results
numFold = 20;
varExplained = 80;
decoderClass = phonemeDecoderClass(numFold,varExplained);
decodeResultStruct = baseDecoder(decoderClass,ieegHGRespNorm,phonemeTrial,'phoneme',[-0.25 0.25]);

 %% Phoneme decoding optimal time-window

decoderObj = decoderClass(10,[10:10:90],1);
timeWin = 0.05:0.05:0.5;
accWin = [];
for iTime = 1:length(timeWin)
    iTime
    for iTer = 1:10
        decodeResultStruct = decoderObj.baseClassify(ieegHGRespNorm,phonemeTrial.phonemeUnit(:,3)',[-timeWin(iTime) timeWin(iTime)], sigChannel);
        accWin(iTime,iTer) = decodeResultStruct.accPhonemeUnBias;
    end
end
%%
w95 = 0.7193;
accWin = squeeze(accWinAll(1,:,:));
 h = figure('Units', 'pixels', ...
    'Position', [100 100 750 500]);
set(h, 'DefaultTextFontSize', 15);

h = boxplot(accWin'.*100', 'symbol','','Whisker',w95,'Colors','k');
hold on;
plot(1:10,median(accWin,2).*100,'LineWidth',2);

%yline(mean(decompStructClassChance.silScore(:)), '--','chance','LineWidth',1);
set(h,{'linew'},{1})
axis square
ylim([0 70])
yline(11.11, '--','chance','LineWidth',1,'Color','r');
set(gca,'XTickLabels', timeWin.*2000);
xlabel('Window size (ms)')

ylabel('P1-Accuracy (%)');

% ylim([0 35])
formatTicks(gca)


%% Phoneme maps
decodeResultStruct = indChanClassify(decoderObj,ieegHGRespNorm,phonemeTrial,'class',[-0.5 0.5]);
chanLabAuc = [];
for iChan = 1:length(decodeResultStruct)
    chanLabAuc(iChan,:) = decodeResultStruct{iChan}.aucAll(1,:); 
end
labels = {'low','high','labials','dorsals'};
%labels = {'a','ae','i','u','b','p','v','g','k'};

%chanLabDiff = abs(chanLabAuc - mean(chanLabAuc,2));
chanLabDiff = (chanLabAuc - 0.5);
figure;

for iPhon = 1:4
     subplot(2,2,iPhon)
%figure
    chanVal = chanView(squeeze(chanLabAuc(:,iPhon)),(chanMap),selectedChannels,[],labels{iPhon},[0.6 1],[],0.25);
    
    colormap(flipud(bone(4096)));
    prc4cutoff = 90;
    cutoff = prctile(chanVal(:),prc4cutoff)
    hold on;
    
    
    chanValMark = true(size(chanVal));
    chanValMark(chanVal<cutoff) = false;
    chanValTemp = chanVal;
    chanValTemp(chanVal<cutoff) = 0;
    measurements = regionprops(chanValMark, chanValTemp, 'all');
    
    
   
    set(gca,'FontSize',15);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis tight
    if(isempty(measurements)||round(cutoff,2)<0.75)
        continue;
    else
     [intensePhoneme(iPhon),maxId] = max([measurements(:).MaxIntensity]);
    hold on;
    scatter(measurements(maxId).ConvexHull(:,1)',measurements(maxId).ConvexHull(:,2)',20,'r','filled')
    x = sqrt(measurements(maxId).ConvexArea);
    phonArea = (1.33*x)^2
    end
    centerOfMass = measurements(maxId).WeightedCentroid;
    scatter(centerOfMass(1),centerOfMass(2),50,'x','LineWidth',2)
    
end

%% Covariance analysis
elecComb = nchoosek(1:length(selectedChannels),2);
eucDist = [];
pitch = 1.33;
sigCov = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = pitch.*sqrt((x2-x1)^2+(y2-y1)^2);
    for iPhoneme = 1:4
        gamma1 = squeeze(chanLabAuc((elecComb(eD,1)),iPhoneme));
        gamma2 = squeeze(chanLabAuc((elecComb(eD,2)),iPhoneme));
        covtemp = cov(gamma1,gamma2);
        sigCov(eD,iPhoneme) = sqrt(-2*gamma1*gamma2+(gamma1^2+gamma2^2));
        %sigCov(eD,iPhoneme) = sqrt(gamma1*gamma2);
        
    end
end

 for iPhoneme = 1:4
    eucDistUnique = unique(eucDist);
    for eDU = 1:length(eucDistUnique)
        eIDs = find(eucDist == eucDistUnique(eDU));
        sigCovMean(eDU,iPhoneme) = median(sigCov(eIDs,iPhoneme));
        sigCovStd(eDU,iPhoneme) = std(sigCov(eIDs,iPhoneme))./sqrt(length(eIDs));   
    end
 end
% figure;
% errorbar(eucDistUnique, sigCovMean(:,2)',sigCovStd(:,2)');
% hold on;
% errorbar(eucDistUnique, sigCovMean(:,7)',sigCovStd(:,7)');
% errorbar(eucDistUnique, sigCovMean(:,11)',sigCovStd(:,11)');
% errorbar(eucDistUnique, sigCovMean(:,12)',sigCovStd(:,12)');

figure;
for iPhoneme = 1:4
errorbar(eucDistUnique, sigCovMean(:,iPhoneme)',sigCovStd(:,iPhoneme)');
hold on;
end

axis square
xlabel('Distance (mm)');
ylabel('Semivariance');
set(gca,'FontSize',20);
xlim([0 10])
legend('low','high','labial','dorsal');

%% Extract Spatial Smooth
accSmooth = zeros(10,8);
for iWindow = 1:8
    iWindow
    if(iWindow == 1)
        ieegHGRespNorm = ieegResponseCar.extractHiGammaNorm(ieegBaseCar,200,[-0.5 0.5],[-0.5 0]);
        
    else
        ieegBaseCarSmooth = ieegBaseCar.spatialSmoothMicro([iWindow iWindow]);
        ieegResponseCarSmooth = ieegResponseCar.spatialSmoothMicro([iWindow iWindow]);
        ieegHGRespNorm = ieegResponseCarSmooth.extractHiGammaNorm(ieegBaseCarSmooth,200,[-0.5 0.5],[-0.5 0]);   
        
    end
    for iTer = 1:10
        decoderClass = phonemeDecoderClass(20,80);
        decodeResultStruct = baseClassify(decoderClass,ieegHGRespNorm,phonemeTrial,'phoneme',[-0.25 0.25]);
        accSmooth(iTer,iWindow) = decodeResultStruct.accPhonemeUnBias(3);
    end
end



elecSpaceSampStr = [];
accSpaceSampAll = [];
phonErrorSpaceSampAll = [];
for iWindow = 1:8  
    elecSpaceSampStrTemp = [num2str(iWindow) 'x' num2str(iWindow)];
    elecSpaceSampStr = [elecSpaceSampStr; repmat({elecSpaceSampStrTemp},10,1)];
    accSpaceSampAll = [accSpaceSampAll accSmooth(:,iWindow)'];
    accSpaceMean(iWindow) = median(accSmooth(:,iWindow)');
   
end 
save('PhonemeDecodeSpaceSmoothThirdPhoneme.mat','accSpaceSampAll','phonErrorSpaceSampAll','elecSpaceSampStr', 'accSpaceMean');

 h = figure('Units', 'pixels', ...
    'Position', [100 100 750 500]);
set(h, 'DefaultTextFontSize', 15); h = boxplot(accSpaceSampAll.*100,elecSpaceSampStr','symbol','','Colors','k');
set(h,{'linew'},{2})
ylim([0 0.7.*100]);
hold on;
plot(1:8,accSpaceMean.*100,'LineWidth',2);
yline(0.1111.*100, '--','chance','LineWidth',2, 'LabelHorizontalAlignment','left');
yline(max(accSpaceMean.*100).*0.95, ':','95% threshold','LineWidth',2, 'Color','r', 'LabelHorizontalAlignment','left');
xlabel('Spatial-average grid size');
ylabel('P1 accuracy (%)');
set(gca,'FontSize',15);
axis square
set(gca,'FontSize',20);


axis square

set( gca                       , ...
    'FontName'   , 'Arial' );


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );

%% HG space-time activation display
ieegGammaPhoneme = [];
labels = {'a','ae','i','u','b','p','v','g','k'};
 phonemeUnit = phonemeTrial.phonemeUnit(:,1)';   
 timeGamma = linspace(ieegHGRespNorm.tw(1),ieegHGRespNorm.tw(2),size(ieegHGRespNorm.data,3));
for iPhoneme = 1:9

    phonIds = (phonemeUnit==iPhoneme);
    ieegGammaMeanTemp = squeeze(mean(ieegHGRespNorm.data(:,phonIds,:),2));
    %travellingWaveMovie(ieegGammaMean,fsD,chanMap,selectedChannelsClean,timeGamma,[-0.25 1],[-1.5 1.5],120,[finger{iFinger} '-clean'],'z-score');

    ieegGammaPhoneme(:,iPhoneme,:) = ieegGammaMeanTemp ;
end


tPlot = [-0.6:0.2:0.6];
phonemeComboLabel =labels;
figure;
tiledlayout(9,length(tPlot),'TileSpacing','compact');


    
    
        for iPhoneme = 1:9
            for iTime = 1:length(tPlot)
        %subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        nexttile;
        val2disp = squeeze(ieegGammaPhoneme(:,iPhoneme,find(timeGamma>=tPlot(iTime),1)));
        chanView(val2disp,chanMap,titl = '',cval = [0 3], isSmooth =0.5);%[1,8,8,16]
        
         
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        colormap((parula(4096)));
        if(iTime==1)
            ylh = ylabel((phonemeComboLabel(iPhoneme)));             
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','top','HorizontalAlignment', 'right')
        end
        if(iPhoneme==1)
            title([num2str(round(tPlot(iTime),2)) ])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',10);
            end
        end