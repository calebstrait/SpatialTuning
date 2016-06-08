function SpatialTuningPaperRasters
close all

figures = [4 5];
%Fig3: offer 1 contra/ipsi
%Fig4: offer 1 contra/ipsi
%Fig5: offer 1 contra/ipsi; value 1 high/low; OFC only

% MAKE RASTERS
for f = figures
    if f ~= 5,
        areas = [1 2 3 4 5];
    else
        areas = 3;
    end
    for a = areas
        if a == 1,
            c = 126;
            load /Users/Jessica/Documents/Writing/GradSchool/Data/StagOps/CESdata.mat CESdata
            spikes = CESdata{c}.spikes;
            sideOfFirst = CESdata{c}.sideOfFirst;
            TS_option1_ON = CESdata{c}.TS_option1_ON;
            choseLefts = zeros(size(sideOfFirst));
            choseLefts((CESdata{c}.choice1stvs2nd == 1 & sideOfFirst == 1) | (CESdata{c}.choice1stvs2nd ~= 1 & sideOfFirst ~= 1)) = 1;
        end
        if a == 2,
            c = 33;
            load /Users/Jessica/Documents/Writing/GradSchool/Data/StagOps/VSdata.mat VSdata
            spikes = VSdata{c}.spikes;
            sideOfFirst = VSdata{c}.sideOfFirst;
            TS_option1_ON = VSdata{c}.TS_option1_ON;
            choseLefts = zeros(size(sideOfFirst));
            choseLefts((VSdata{c}.choice1stvs2nd == 1 & sideOfFirst == 1) | (VSdata{c}.choice1stvs2nd ~= 1 & sideOfFirst ~= 1)) = 1;
        end
        if a == 3,
            if f == 3 || f == 4
                c = 56;
                OFCraw = load('/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/OFC_raw/B130313.1S.plx.mat');
                spikes = OFCraw.SPK01b;
            else
                c = 31;
                OFCraw = load('/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/OFC_raw/B130220.2S.plx.mat');
                spikes = OFCraw.SPK01a;
            end
            load('/Users/Jessica/Documents/Writing/GradSchool/Data/Curiosity/ofc_data.mat');
            OFCdata = data;
            clear data
            sideOfFirst = OFCdata(c).var(9,:)';
            strobes = getstrobes(OFCraw.EVT01,OFCraw.EVT02,OFCraw.EVT03,OFCraw.EVT04,OFCraw.EVT05,OFCraw.EVT06);
            strobeTs = strobes(:,1);
            strobeVs = strobes(:,2);
            TS_option1_ON = strobeTs(strobeVs == 5251);
            choseLefts = OFCdata(c).var(1,:)';
            evL = OFCdata(c).var(5,:)';
            evR = OFCdata(c).var(4,:)';
            offer1Value = zeros(size(evL));
            h = OFCdata(c).var(9,:)';
            offer1Value(h == 1) = evL(h == 1);
            offer1Value(h ~= 1) = evR(h ~= 1);
        end
        if a == 4,
            c = 88;
            loaded = load('/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/dACC_data.mat');
            dACCdata = loaded.data;
            [~,name,~] = fileparts(dACCdata{c}.file);
            dACCraw = load(['/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/dACC_raw/' name '.mat']);
            eval(['spikes = dACCraw.' dACCdata{c}.cell ';']);
            o = dACCdata{c}.vars(:,15);
            h = ones(size(o));
            h(o == 2) = 0;
            sideOfFirst = h;
            strobeTs = dACCraw.Strobed(:,1);
            strobeVs = dACCraw.Strobed(:,2);
            TS_option1_ON = strobeTs(strobeVs == 6020);
            chos1 = dACCdata{c}.vars(:,16);
            chosL = zeros(size(chos1));
            chosL(chos1 == 1 & f == 1 | chos1 ~= 1 & f ~= 1) = 1;
            choseLefts = chosL;
        end
        if a == 5,
            if f == 3
                c = 33;
            else
                c = 67;
            end
            loaded = load('/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/sgACC_data.mat');
            sgACCdata = loaded.data;
            [~,name,~] = fileparts(sgACCdata{c}.file);
            if f == 3
                sgACCraw = load('/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/sgACC_raw/B140902.PTR25.1S.pl2.mat');
            else
                cd /Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/sgACC_raw/
                sgACCraw = load([ name '.mat']);
            end
            eval(['spikes = sgACCraw.' sgACCdata{c}.cell ';']);
            o = sgACCdata{c}.vars(:,15);
            h = ones(size(o));
            h(o == 2) = 0;
            sideOfFirst = h;
            strobeTs = sgACCraw.Strobed(:,1);
            strobeVs = sgACCraw.Strobed(:,2);
            TS_option1_ON = strobeTs(strobeVs == 6020);
            chos1 = sgACCdata{c}.vars(:,16);
            chosL = zeros(size(chos1));
            chosL(chos1 == 1 & f == 1 | chos1 ~= 1 & f ~= 1) = 1;
            choseLefts = chosL;
        end
        
        clearvars -except spikes sideOfFirst TS_option1_ON a c areas f choseLefts offer1Value
        
        ntrials = length(sideOfFirst);
        firstOnLefts = sideOfFirst;
        nSet = [0 0 0 0];
        count = [0 0 0 0];
        setsX = zeros(ntrials,1);
        fig = figure;
        hold on
        for jj = 1:ntrials %count in each set
            if f == 3
                if firstOnLefts(jj) == 1
                    setsX(jj) = 1;
                    nSet(1) = nSet(1) + 1;
                else
                    setsX(jj) = 2;
                    nSet(2) = nSet(2) + 1;
                end
            elseif f == 4
                if choseLefts(jj) == 1
                    setsX(jj) = 1;
                    nSet(1) = nSet(1) + 1;
                else
                    setsX(jj) = 2;
                    nSet(2) = nSet(2) + 1;
                end
            else
                if firstOnLefts(jj) == 1 && offer1Value(jj) >= median(offer1Value)
                    setsX(jj) = 1;
                    nSet(1) = nSet(1) + 1;
                elseif firstOnLefts(jj) == 1 && offer1Value(jj) < median(offer1Value)
                    setsX(jj) = 2;
                    nSet(2) = nSet(2) + 1;
                elseif firstOnLefts(jj) ~= 1 && offer1Value(jj) >= median(offer1Value)
                    setsX(jj) = 3;
                    nSet(3) = nSet(3) + 1;
                elseif firstOnLefts(jj) ~= 1 && offer1Value(jj) < median(offer1Value)
                    setsX(jj) = 4;
                    nSet(4) = nSet(4) + 1;
                end
            end
        end
        for jj = 1:ntrials
            
            [rel_spikes1,rel_spikes2] = deal([]);

            starttime = TS_option1_ON - 1;
            stoptime = TS_option1_ON + 4;
            [rel_spikes3,rel_spikes4] = deal([]);
            if setsX(jj) == 1
                rel_spikes1 = spikes((spikes > starttime(jj)) & (spikes < stoptime(jj)));
                rel_spikes1 = rel_spikes1 - starttime(jj);
                count(1) = count(1) + 1;
            elseif setsX(jj) == 2                
                rel_spikes2 = spikes((spikes > starttime(jj)) & (spikes < stoptime(jj)));
                rel_spikes2 = rel_spikes2 - starttime(jj);
                count(2) = count(2) + 1;
            elseif setsX(jj) == 3
                rel_spikes3 = spikes((spikes > starttime(jj)) & (spikes < stoptime(jj)));
                rel_spikes3 = rel_spikes3 - starttime(jj);
                count(3) = count(3) + 1;
            else
                rel_spikes4 = spikes((spikes > starttime(jj)) & (spikes < stoptime(jj)));
                rel_spikes4 = rel_spikes4 - starttime(jj);
                count(4) = count(4) + 1;
            end
            
            for ii = 1:length(rel_spikes1) % for every spike
                time = rel_spikes1(ii);
                line([time time+.001],[count(1)-0.5 count(1)+0.5],'Color',[0 0 .75]);
                % draw a black vertical line of length 1 at time t (x) and at trial jj (y)
            end
            
            for ii = 1:length(rel_spikes2)
                time = rel_spikes2(ii);
                line([time time+.001],[nSet(1)+count(2)-0.5 nSet(1)+count(2)+0.5],'Color',[0 .75 0]);
            end
            
            for ii = 1:length(rel_spikes3)
                time = rel_spikes3(ii);
                line([time time+.001],[nSet(1)+nSet(2)+count(3)-0.5 nSet(1)+nSet(2)+count(3)+0.5],'Color',[.75 .75 0]);
            end
            
            for ii = 1:length(rel_spikes4)
                time = rel_spikes4(ii);
                line([time time+.001],[nSet(1)+nSet(2)+nSet(3)+count(4)-0.5 nSet(1)+nSet(2)+nSet(3)+count(4)+0.5],'Color',[.75 0 .75]);
            end
            
        end
        names = {'VM','VS','OF','DA','SA'};
        vlineCords = [1 2];
        title(sprintf('%s Cell %i',names{a},c),'FontSize',14);
        xlabel('Time (ms)'); % Time is in millisecond
        ylabel('Trial number');
        vlineColor(vlineCords(1),[.5 .5 .5]);
        vlineColor(vlineCords(2),[.5 .5 .5]);
        axis([0 5 0 ntrials]);
        
        hold off
        saveas(fig, ['/Users/Jessica/Documents/Writing/GradSchool/Papers/JNeurophysiology/Raw_Rasters/' num2str(f) names{a} '.eps'],'eps2c');
        close all
    end
end

cd /Users/Jessica/Desktop/
keyboard
end