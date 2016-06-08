function SpatialTuningPaperGLM

%%%%%%LOAD DATA
load /Users/cstrait/Documents/Writing/GradSchool/Data/StagOps/vmPFCdata.mat
load /Users/cstrait/Documents/Writing/GradSchool/Data/StagOps/VSgamblingdata.mat
VSdata = VSgamblingdata;
clear VSgamblingdata
load /Users/cstrait/Desktop/OFC2015/dACC_data.mat
dACCdata = data;
clear data
load /Users/cstrait/Desktop/OFC2015/sgACC_data.mat
sgACCdata = data;
clear data
load /Users/cstrait/Documents/Writing/GradSchool/Data/Curiosity/ofc_data.mat
OFCdata = data;
clear data

%%%%%FORMAT DATA
[VM_E1,VM_E2,VM_E3,VS_E1,VS_E2,VS_E3,DA_E1,DA_E2,DA_E3,SA_E1,SA_E2,SA_E3,OF_E1,OF_E2,OF_E3,CellNum,FirstOnLeft,ChoseLeft,ChosenRewardSize] = deal([]);
trialCounts = zeros(5, 1);
for i = 1:length(vmPFCdata)
    VM_E1 = [VM_E1; mean(vmPFCdata{i}.psth(:,150:174),2)];
    VM_E2 = [VM_E2; mean(vmPFCdata{i}.psth(:,200:224),2)];
    VM_E3 = [VM_E3; mean(vmPFCdata{i}.psth(:,265:289),2)];
    if isempty(CellNum)
        CellNum = ones(length(vmPFCdata{i}.vars)*3,1);
    else
        CellNum = [CellNum; (max(CellNum)+1) * ones(length(vmPFCdata{i}.vars)*3,1)];
    end
    FirstOnLeft = [FirstOnLeft; vmPFCdata{i}.vars(:,7)];
    ChoseLeft = [ChoseLeft; vmPFCdata{i}.vars(:,8)];
    cho1vs2 = vmPFCdata{i}.vars(:,9);
    ev1 = vmPFCdata{i}.vars(:,3);
    ev2 = vmPFCdata{i}.vars(:,6);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    ChosenRewardSize = [ChosenRewardSize; choRewSize];
    trialCounts(1) = trialCounts(1) + length(vmPFCdata{i}.vars);
end
for i = 1:length(VSdata)
    VS_E1 = [VS_E1; mean(VSdata{i}.psth(:,150:174),2)];
    VS_E2 = [VS_E2; mean(VSdata{i}.psth(:,200:224),2)];
    VS_E3 = [VS_E3; mean(VSdata{i}.psth(:,265:289),2)];
    CellNum = [CellNum; (max(CellNum)+1) * ones(length(VSdata{i}.vars)*3,1)];
    FirstOnLeft = [FirstOnLeft; VSdata{i}.vars(:,7)];
    ChoseLeft = [ChoseLeft; VSdata{i}.vars(:,8)];
    cho1vs2 = VSdata{i}.vars(:,9);
    ev1 = VSdata{i}.vars(:,3);
    ev2 = VSdata{i}.vars(:,6);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    ChosenRewardSize = [ChosenRewardSize; choRewSize];
    trialCounts(2) = trialCounts(2) + length(VSdata{i}.vars);
end
for i = 1:length(dACCdata)
    DA_E1 = [DA_E1; mean(dACCdata{i}.psth(:,250:274),2)];
    DA_E2 = [DA_E2; mean(dACCdata{i}.psth(:,287:311),2)];
    DA_E3 = [DA_E3; mean(dACCdata{i}.psth(:,340:364),2)];
    CellNum = [CellNum; (max(CellNum)+1) * ones(length(dACCdata{i}.vars)*3,1)];
    o = dACCdata{i}.vars(:,15);
    f = ones(size(o));
    f(o == 2) = 0;
    FirstOnLeft = [FirstOnLeft; f];
    chos1 = dACCdata{i}.vars(:,16);
    chosL = zeros(size(chos1));
    chosL(chos1 == 1 & f == 1 | chos1 ~= 1 & f ~= 1) = 1;
    ChoseLeft = [ChoseLeft; chosL];
    clear o chos1
    cho1vs2 = dACCdata{i}.vars(:,16);
    evL = dACCdata{i}.vars(:,5) .* dACCdata{i}.vars(:,12) + dACCdata{i}.vars(:,6) .* (1-dACCdata{i}.vars(:,12));
    evR = dACCdata{i}.vars(:,7) .* dACCdata{i}.vars(:,13) + dACCdata{i}.vars(:,8) .* (1-dACCdata{i}.vars(:,13));
    [ev1,ev2] = deal(zeros(size(evL)));
    ev1(f == 1) = evL(f == 1);
    ev1(f ~= 1) = evR(f ~= 1);
    ev2(f == 1) = evR(f == 1);
    ev2(f ~= 1) = evL(f ~= 1);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    ChosenRewardSize = [ChosenRewardSize; choRewSize];
    trialCounts(3) = trialCounts(3) + length(dACCdata{i}.vars);
end
for i = 1:length(sgACCdata)
    SA_E1 = [SA_E1; mean(sgACCdata{i}.psth(:,250:274),2)];
    SA_E2 = [SA_E2; mean(sgACCdata{i}.psth(:,287:311),2)];
    SA_E3 = [SA_E3; mean(sgACCdata{i}.psth(:,340:364),2)];
    CellNum = [CellNum; (max(CellNum)+1) * ones(length(sgACCdata{i}.vars)*3,1)];
    o = sgACCdata{i}.vars(:,15);
    f = ones(size(o));
    f(o == 2) = 0;
    FirstOnLeft = [FirstOnLeft; f];
    chos1 = sgACCdata{i}.vars(:,16);
    chosL = zeros(size(chos1));
    chosL(chos1 == 1 & f == 1 | chos1 ~= 1 & f ~= 1) = 1;
    ChoseLeft = [ChoseLeft; chosL];
    clear o chos1
    cho1vs2 = sgACCdata{i}.vars(:,16);
    evL = sgACCdata{i}.vars(:,5) .* sgACCdata{i}.vars(:,12) + sgACCdata{i}.vars(:,6) .* (1-sgACCdata{i}.vars(:,12));
    evR = sgACCdata{i}.vars(:,7) .* sgACCdata{i}.vars(:,13) + sgACCdata{i}.vars(:,8) .* (1-sgACCdata{i}.vars(:,13));
    [ev1,ev2] = deal(zeros(size(evL)));
    ev1(f == 1) = evL(f == 1);
    ev1(f ~= 1) = evR(f ~= 1);
    ev2(f == 1) = evR(f == 1);
    ev2(f ~= 1) = evL(f ~= 1);
    choRewSize = zeros(size(cho1vs2));
    choRewSize(cho1vs2 == 1) = ev1(cho1vs2 == 1);
    choRewSize(cho1vs2 == 2) = ev2(cho1vs2 == 2);
    ChosenRewardSize = [ChosenRewardSize; choRewSize];
    trialCounts(4) = trialCounts(4) + length(sgACCdata{i}.vars);
end
for i = 1:length(OFCdata)
    OF_E1 = [OF_E1; mean(OFCdata(i).epoch1(:,1:24),2)];
    OF_E2 = [OF_E2; mean(OFCdata(i).epoch2(:,1:24),2)];
    OF_E3 = [OF_E3; mean(OFCdata(i).epoch4(:,1:24),2)];
    CellNum = [CellNum; (max(CellNum)+1) * ones(length(OFCdata(i).var)*3,1)];
    FirstOnLeft = [FirstOnLeft; OFCdata(i).var(9,:)'];
    ChoseLeft = [ChoseLeft; OFCdata(i).var(1,:)'];
    choL = OFCdata(i).var(1,:)';
    choRewSize = zeros(size(choL));
    evL = OFCdata(i).var(5,:)';
    evR = OFCdata(i).var(4,:)';
    choRewSize(choL == 1) = evL(choL == 1);
    choRewSize(choL == 0) = evR(choL == 0);
    ChosenRewardSize = [ChosenRewardSize; choRewSize];
    trialCounts(5) = trialCounts(5) + length(OFCdata(i).var);
end
Spikes = [VM_E1;VM_E2;VM_E3;VS_E1;VS_E2;VS_E3;DA_E1;DA_E2;DA_E3;SA_E1;SA_E2;SA_E3;OF_E1;OF_E2;OF_E3];
Experiment = [[ones(size(VM_E1));ones(size(VM_E2));ones(size(VM_E3))];...
    [ones(size(VS_E1));ones(size(VS_E2));ones(size(VS_E3))];...
    2 * [ones(size(DA_E1));ones(size(DA_E2));ones(size(DA_E3))];...
    2 * [ones(size(SA_E1));ones(size(SA_E2));ones(size(SA_E3))];...
    3 * [ones(size(OF_E1));ones(size(OF_E2));ones(size(OF_E3))]];
Region = [[ones(size(VM_E1));ones(size(VM_E2));ones(size(VM_E3))];...
    2 * [ones(size(VS_E1));ones(size(VS_E2));ones(size(VS_E3))];...
    3 * [ones(size(DA_E1));ones(size(DA_E2));ones(size(DA_E3))];...
    4 * [ones(size(SA_E1));ones(size(SA_E2));ones(size(SA_E3))];...
    5 * [ones(size(OF_E1));ones(size(OF_E2));ones(size(OF_E3))]];
Epoch = [ones(size(VM_E1));2*ones(size(VM_E2));3*ones(size(VM_E3));...
    ones(size(VS_E1));2*ones(size(VS_E2));3*ones(size(VS_E3));...
    ones(size(DA_E1));2*ones(size(DA_E2));3*ones(size(DA_E3));...
    ones(size(SA_E1));2*ones(size(SA_E2));3*ones(size(SA_E3));...
    ones(size(OF_E1));2*ones(size(OF_E2));3*ones(size(OF_E3))];
FirstOnLeft = [FirstOnLeft(1:trialCounts(1));FirstOnLeft(1:trialCounts(1));FirstOnLeft(1:trialCounts(1));...
    FirstOnLeft(1:trialCounts(2));FirstOnLeft(1:trialCounts(2));FirstOnLeft(1:trialCounts(2));...
    FirstOnLeft(1:trialCounts(3));FirstOnLeft(1:trialCounts(3));FirstOnLeft(1:trialCounts(3));...
    FirstOnLeft(1:trialCounts(4));FirstOnLeft(1:trialCounts(4));FirstOnLeft(1:trialCounts(4));...
    FirstOnLeft(1:trialCounts(5));FirstOnLeft(1:trialCounts(5));FirstOnLeft(1:trialCounts(5))];
ChoseLeft = [ChoseLeft(1:trialCounts(1));ChoseLeft(1:trialCounts(1));ChoseLeft(1:trialCounts(1));...
    ChoseLeft(1:trialCounts(2));ChoseLeft(1:trialCounts(2));ChoseLeft(1:trialCounts(2));...
    ChoseLeft(1:trialCounts(3));ChoseLeft(1:trialCounts(3));ChoseLeft(1:trialCounts(3));...
    ChoseLeft(1:trialCounts(4));ChoseLeft(1:trialCounts(4));ChoseLeft(1:trialCounts(4));...
    ChoseLeft(1:trialCounts(5));ChoseLeft(1:trialCounts(5));ChoseLeft(1:trialCounts(5))];
ChosenRewardSize = [ChosenRewardSize(1:trialCounts(1));ChosenRewardSize(1:trialCounts(1));ChosenRewardSize(1:trialCounts(1));...
    ChosenRewardSize(1:trialCounts(2));ChosenRewardSize(1:trialCounts(2));ChosenRewardSize(1:trialCounts(2));...
    ChosenRewardSize(1:trialCounts(3));ChosenRewardSize(1:trialCounts(3));ChosenRewardSize(1:trialCounts(3));...
    ChosenRewardSize(1:trialCounts(4));ChosenRewardSize(1:trialCounts(4));ChosenRewardSize(1:trialCounts(4));...
    ChosenRewardSize(1:trialCounts(5));ChosenRewardSize(1:trialCounts(5));ChosenRewardSize(1:trialCounts(5))];
clearvars -except Spikes CellNum Experiment Region FirstOnLeft ChoseLeft ChosenRewardSize Epoch

%%%%%GLM

%%% Reviewer:
% offer position, choice position, and reward value. 

%%% Mathworks:
% [pvals,tbl,stats] = anovan(mileage, {factory carmod}, ...
% 'model',2, 'random',1,'varnames',{'Factory' 'Car Model'});

disp('ANOVAing');
[pvals,tbl,stats] = anovan(Spikes, {CellNum Experiment Region FirstOnLeft ChoseLeft ChosenRewardSize Epoch}, ...
    'random',[1 2],'varnames',{'CellNum' 'Experiment' 'Region' 'FirstOnLeft' 'ChoseLeft' 'ChosenRewardSize' 'Epoch'});
disp('Done!');

% % X = [ones(406,1) Acceleration Horsepower];
% % Z = {ones(406,1),Acceleration};
% % G = {Model_Year,Model_Year};
% % lme = fitlmematrix(X,Y,Z,G,'FixedEffectPredictors',....
% % {'Intercept','Acceleration','Horsepower'},'RandomEffectPredictors',...
% % {{'Intercept','Acceleration'}},'RandomEffectGroups',{'Model_Year'});
% 
% % Fixed effects matrix
% X = [ones(length(FirstOnLeft),1) Region FirstOnLeft ChoseLeft ChosenRewardSize Epoch];
% 
% % Random effects matrix
% Z = {ones(length(FirstOnLeft),1), CellNum, Experiment};
% 
% lme = fitlmematrix(X,Spikes,Z,[],...
%     'FixedEffectPredictors',{'Intercept','Region','FirstOnLeft','ChoseLeft','ChosenRewardSize','Epoch'},...
%     'RandomEffectPredictors',{{'Intercept'},{'CellNum'},{'Experiment'}})

keyboard
end