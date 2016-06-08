function SpatialTuningEyeAnalysis
close all

%%%%%%PREP STRUCTURE
C = {5,1};
for i = 1:5
    C{i}.FirstOnLeft = {};
    C{i}.ChoseLeft = {};
    C{i}.continuous = {};
end

%%%%%%LOAD DATA
load /Users/Jessica/Desktop/Writing/GradSchool/Data/StagOps/vmPFCdata.mat
load /Users/Jessica/Desktop/Writing/GradSchool/Data/StagOps/VSgamblingdata.mat
VSdata = VSgamblingdata;
clear VSgamblingdata
load /Users/Jessica/Desktop/JNeurophysiology/dACC_data.mat
dACCdata = data;
clear data
load /Users/Jessica/Desktop/JNeurophysiology/sgACC_data.mat
sgACCdata = data;
clear data
load /Users/Jessica/Desktop/Writing/GradSchool/Data/Curiosity/ofc_data.mat
OFCdata = data;
clear data

%%%%%%LOAD EYE DATA
for c=1:length(dACCdata) %dACC
    [~,NAME,EXT] = fileparts(dACCdata{c}.file);
    C{3}.HaveEyeData{c} = 0;
    if exist(['/Users/Jessica/Desktop/JNeurophysiology/Eye_data/dACC/' NAME(1:7) '.Eye.PTACC' NAME(8:end) EXT], 'file') ~=0
        theseData = load(['/Users/Jessica/Desktop/JNeurophysiology/Eye_data/dACC/' NAME(1:7) '.Eye.PTACC' NAME(8:end) EXT]);
        eval(['continuous.spikes = theseData.' dACCdata{c}.cell ';']);
        continuous.eyeX = theseData.FP13;
        continuous.eyeY = theseData.FP14;
        continuous.eyeX_TS = theseData.FP13_ts(end);
        continuous.eyeY_TS = theseData.FP14_ts(end);
        C{3}.continuous{c} = continuous;
        theseData.Strobed = trimStrobedACC(theseData.Strobed);
        C{3}.Op1On_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 6030),1);
        C{3}.Op1Of_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 6040),1);
        C{3}.BothOn_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 6080),1);
        found = find(theseData.Strobed(:,2) == 6010);
        C{3}.FirstOnLeft{c} = theseData.Strobed(found+14,2);
        C{3}.ChoseLeft{c} = theseData.Strobed(found+15,2);
        if ~(length(C{3}.Op1On_TS{c}) ~= length(C{3}.BothOn_TS{c}) || length(C{3}.Op1On_TS{c}) ~= length(C{3}.ChoseLeft{c}) || length(C{3}.FirstOnLeft{c}) ~= length(C{3}.ChoseLeft{c}))
            C{3}.HaveEyeData{c} = 1;
        end
    end
end
for c=1:length(sgACCdata) %sgACC
    [~,NAME,EXT] = fileparts(sgACCdata{c}.file);
    C{4}.HaveEyeData{c} = 0;
    if exist(['/Users/Jessica/Desktop/JNeurophysiology/Eye_data/sgACC/' NAME(1:7) '.Eye.PTR25.' NAME(15:16) EXT], 'file') == 2
        theseData = load(['/Users/Jessica/Desktop/JNeurophysiology/Eye_data/sgACC/' NAME(1:7) '.Eye.PTR25.' NAME(15:16) EXT]);
        eval(['continuous.spikes = theseData.' sgACCdata{c}.cell ';']);
        continuous.eyeX = theseData.EyeX;
        continuous.eyeY = theseData.EyeY;
        continuous.eyeX_TS = theseData.FP13_ts(end);
        continuous.eyeY_TS = theseData.FP14_ts(end);
        C{4}.continuous{c} = continuous;
        theseData.Strobed = trimStrobedACC(theseData.Strobed);
        C{4}.Op1On_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 6030),1);
        C{4}.Op1Of_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 6040),1);
        C{4}.BothOn_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 6080),1);
        found = find(theseData.Strobed(:,2) == 6010);
        C{4}.FirstOnLeft{c} = theseData.Strobed(found+14,2);
        C{4}.ChoseLeft{c} = theseData.Strobed(found+15,2);
        if ~(length(C{4}.Op1On_TS{c}) ~= length(C{4}.BothOn_TS{c}) || length(C{4}.Op1On_TS{c}) ~= length(C{4}.ChoseLeft{c}) || length(C{4}.FirstOnLeft{c}) ~= length(C{4}.ChoseLeft{c}))
            C{4}.HaveEyeData{c} = 1;
        end
    end
end
for c=1:length(OFCdata) %OFC            Must replace Joker data with Hobbes data and check for more Batman data
    [filename,letter] = ACUR_CellInfo(c);
    C{5}.HaveEyeData{c} = 0;
    if exist(['/Users/Jessica/Desktop/JNeurophysiology/Eye_data/OFC/' filename], 'file') == 2
        theseData = load(['/Users/Jessica/Desktop/JNeurophysiology/Eye_data/OFC/' filename]);
        eval(['continuous.spikes = theseData.SPK0' letter ';']);
        continuous.eyeX = theseData.Eye_X;
        continuous.eyeY = theseData.Eye_Y;
        continuous.eyeX_TS = theseData.FP13_ts(end);
        continuous.eyeY_TS = theseData.FP14_ts(end);
        C{5}.continuous{c} = continuous;
        theseData.Strobed = getstrobes(theseData.EVT01,theseData.EVT02,theseData.EVT03,theseData.EVT04,theseData.EVT05,theseData.EVT06);
        theseData.Strobed = trimStrobedOFC(theseData.Strobed);
        C{5}.Op1On_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 5251),1);
        C{5}.Op1Of_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 5252),1);
        C{5}.BothOn_TS{c} = theseData.Strobed((theseData.Strobed(:,2) == 5300),1);
        FoL = theseData.Strobed(theseData.Strobed(:,2) == 5255 | theseData.Strobed(:,2) == 5256,2);
        FoL(FoL == 5255) = 1; FoL(FoL == 5256) = 0;
        C{5}.FirstOnLeft{c} = FoL;
        CL = theseData.Strobed(theseData.Strobed(:,2) == 5060 | theseData.Strobed(:,2) == 5070,2);
        CL(CL == 5060) = 1; CL(CL == 5070) = 0;
        C{5}.ChoseLeft{c} = CL;
        C{5}.HaveEyeData{c} = 1;
    end
end

%%%%%FORMAT DATA
for i = 1:length(vmPFCdata)
end
for i = 1:length(VSdata)
end
for i = 1:length(dACCdata)
    if C{3}.HaveEyeData{i} == 1
        cont = C{3}.continuous{i};
        [EP1_EyeX,EP1_EyeY,EP3_EyeX,EP3_EyeY,EP1_spikes,EP3_spikes] = deal(zeros(size(C{3}.BothOn_TS{i})));
        for j = 1:length(C{3}.Op1On_TS{i})
            startEp1 = C{3}.Op1On_TS{i};
            endEp1 = C{3}.Op1Of_TS{i};
            startEp3 = C{3}.BothOn_TS{i};
            endEp3 = startEp3 + 1;
            endIndex = min(round(endEp3(j) / 0.001),length(cont.eyeX));
            EP1_EyeX(j) = mean(cont.eyeX(round(startEp1(j) / 0.001):round(endEp1(j) / 0.001)));
            EP1_EyeY(j) = mean(cont.eyeY(round(startEp1(j) / 0.001):round(endEp1(j) / 0.001)));
            EP3_EyeX(j) = mean(cont.eyeX(round(startEp3(j) / 0.001):endIndex));
            EP3_EyeY(j) = mean(cont.eyeY(round(startEp3(j) / 0.001):endIndex));
            EP1_spikes(j) = sum(cont.spikes(cont.spikes>=startEp1(j) & cont.spikes<endEp1(j)));
            EP3_spikes(j) = sum(cont.spikes(cont.spikes>=startEp3(j) & cont.spikes<endEp3(j)));
        end
        C{3}.EP1_EyeX{i} = EP1_EyeX;
        C{3}.EP1_EyeY{i} = EP1_EyeY;
        C{3}.EP3_EyeX{i} = EP3_EyeX;
        C{3}.EP3_EyeY{i} = EP3_EyeY;
        C{3}.EP1_spikes{i} = EP1_spikes;
        C{3}.EP3_spikes{i} = EP3_spikes;
    end
end
for i = 1:length(sgACCdata)
    if C{4}.HaveEyeData{i} == 1
        cont = C{4}.continuous{i};
        [EP1_EyeX,EP1_EyeY,EP3_EyeX,EP3_EyeY,EP1_spikes,EP3_spikes] = deal(zeros(size(C{4}.BothOn_TS{i})));
        for j = 1:length(C{4}.Op1On_TS{i})
            startEp1 = C{4}.Op1On_TS{i};
            endEp1 = C{4}.Op1Of_TS{i};
            startEp3 = C{4}.BothOn_TS{i};
            endEp3 = startEp3 + 1;
            endIndex = min(round(endEp3(j) / 0.001),length(cont.eyeX));
            EP1_EyeX(j) = mean(cont.eyeX(round(startEp1(j) / 0.001):round(endEp1(j) / 0.001)));
            EP1_EyeY(j) = mean(cont.eyeY(round(startEp1(j) / 0.001):round(endEp1(j) / 0.001)));
            EP3_EyeX(j) = mean(cont.eyeX(round(startEp3(j) / 0.001):endIndex));
            EP3_EyeY(j) = mean(cont.eyeY(round(startEp3(j) / 0.001):endIndex));
            EP1_spikes(j) = sum(cont.spikes(cont.spikes>=startEp1(j) & cont.spikes<endEp1(j)));
            EP3_spikes(j) = sum(cont.spikes(cont.spikes>=startEp3(j) & cont.spikes<endEp3(j)));
        end
        C{4}.EP1_EyeX{i} = EP1_EyeX;
        C{4}.EP1_EyeY{i} = EP1_EyeY;
        C{4}.EP3_EyeX{i} = EP3_EyeX;
        C{4}.EP3_EyeY{i} = EP3_EyeY;
        C{4}.EP1_spikes{i} = EP1_spikes;
        C{4}.EP3_spikes{i} = EP3_spikes;
    end
end
for i = 1:length(OFCdata)
    if C{5}.HaveEyeData{i} == 1
        cont = C{5}.continuous{i};
        [EP1_EyeX,EP1_EyeY,EP3_EyeX,EP3_EyeY,EP1_spikes,EP3_spikes] = deal(zeros(size(C{5}.BothOn_TS{i})));
        for j = 1:length(C{5}.Op1On_TS{i})
            startEp1 = C{5}.Op1On_TS{i};
            endEp1 = C{5}.Op1Of_TS{i};
            startEp3 = C{5}.BothOn_TS{i};
            endEp3 = startEp3 + 1;
            endIndex = min(round(endEp3(j) / 0.001),length(cont.eyeX));
            EP1_EyeX(j) = mean(cont.eyeX(round(startEp1(j) / 0.001):round(endEp1(j) / 0.001)));
            EP1_EyeY(j) = mean(cont.eyeY(round(startEp1(j) / 0.001):round(endEp1(j) / 0.001)));
            EP3_EyeX(j) = mean(cont.eyeX(round(startEp3(j) / 0.001):endIndex));
            EP3_EyeY(j) = mean(cont.eyeY(round(startEp3(j) / 0.001):endIndex));
            EP1_spikes(j) = sum(cont.spikes(cont.spikes>=startEp1(j) & cont.spikes<endEp1(j)));
            EP3_spikes(j) = sum(cont.spikes(cont.spikes>=startEp3(j) & cont.spikes<endEp3(j)));
        end
        C{5}.EP1_EyeX{i} = EP1_EyeX;
        C{5}.EP1_EyeY{i} = EP1_EyeY;
        C{5}.EP3_EyeX{i} = EP3_EyeX;
        C{5}.EP3_EyeY{i} = EP3_EyeY;
        C{5}.EP1_spikes{i} = EP1_spikes;
        C{5}.EP3_spikes{i} = EP3_spikes;

    end
end

C{1}.leng = length(vmPFCdata);
C{2}.leng = length(VSdata);
C{3}.leng = length(dACCdata);
C{4}.leng = length(sgACCdata);
C{5}.leng = length(OFCdata);

clearvars -except C

%%%%%ANALYSIS
names = {'VM','VS','DA','SA','OF'};
smo = 20;
width = 750;
vlineCords = [250 287 330] - 188; %offer 1, offer 2, both offers on
for a = 3:5 %*%*%*%
    [countX,countY] = deal(zeros(width,1));
    
    dataCount = 0;
    for c = 1:C{a}.leng % # cells
        if C{a}.HaveEyeData{c} == 1
            dataCount = dataCount + 1;
            for x = 1:width
                
                % Are spikes correlated with EyeXposition and EyeYposition?
                C{a}.continuous{c}.eyeX(x:x+24); %*%*%*%*%*%*%*%*
                input = [ones(length(C{a}.EP1_spikes{c}),1) ];
                [~,~,stats] = glmfit(input,C{a}.EP1_spikes{c});
                if stats.p(2) < .05, countX(x) = countX(x) + 1;end
                if stats.p(3) < .05, countY(x) = countY(x) + 1;end
                
            end
        end
    end
    figure
    hold on
    toplotX = smooth(countX / dataCount,smo);
    toplotY = smooth(countY / dataCount,smo);
    plot(toplotX,'b');
    plot(toplotY,'g');
    vlineColor(vlineCords(1),[.5 .5 .5]);
    vlineColor(vlineCords(2),[.5 .5 .5]);
    vlineColor(vlineCords(3),[.5 .5 .5]); %1.6s between offer 1 on and both offers on
    ylabel('%cells','FontSize',14);
    xlabel('time','FontSize',14);
    title(sprintf('%% cells with spikes correlated to eye position %s n=%i/%i',names{a},dataCount,C{a}.leng),'FontSize',14);
    legend({'Eye X','Eye Y'});
    hold off
end

cd /Users/Jessica/Desktop/
keyboard
end

function s = trimStrobedACC(s)
op1onTS = 6030;
op1offTS = 6040;
endTS = 6080;
ITIts = 6010;

starts = find(s(:,2) == op1onTS);
mids = find(s(:,2) == op1offTS);
ends = find(s(:,2) == endTS);
ITIs = find(s(:,2) == ITIts);

for x = [mids(1) ends(1) ITIs(1)]
    if x < starts(1), s(x,2) = -1;end
end
for x = [starts(end) mids(end) ends(end) ITIs(end)]
    if x >= starts(end), s(x,2) = -1;end
end

s((s(:,2) ~= op1onTS & s(:,2) ~= op1offTS & s(:,2) ~= endTS & s(:,2) ~= ITIts),2) = 0; %%%

end

function s = trimStrobedOFC(s)
op1onTS=5251;
op1offTS=5252;
endTS=5300;
choseL=5060;
choseR=5070;
firstOnL=5255;
firstOnR=5256;

starts = find(s(:,2) == op1onTS);
mids = find(s(:,2) == op1offTS);
ends = find(s(:,2) == endTS);
choices = find(s(:,2) == choseL | s(:,2) == choseR);
firstLR = find(s(:,2) == firstOnL | s(:,2) == firstOnR);

for x = [mids(1) ends(1) choices(1) firstLR(1)]
    if x < starts(1), s(x,2) = -1;end
end
for x = [starts(end) mids(end) ends(end) firstLR(end)]
    if x > choices(end), s(x,2) = -1;end
end

end