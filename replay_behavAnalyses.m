%% REPLAY behavioural analyses

clc;clear

% paths
paths           = [];
paths.parent    = '/Users/alex/Documents/replay_dat/';
paths.mri       = [paths.parent 'mri/'];
paths.behav     = [paths.parent 'behav/'];

% experimental info
IDs = {'102'; '103'; '104'; '105'; '106'; '107'; '108'; '109'; '110'; '111'; '112';...
    '114'; '115'; '119'; '120'; '122'};

behav=[]; 
for id=1:length(IDs)
    
    % load
    for enc=1:2
        behav.enc{id,enc}  = load([paths.behav IDs{id} '_enc' num2str(enc) '.mat']);
    end
    for ret1=1:4
        behav.ret1{id,ret1}  = load([paths.behav IDs{id} '_ret1' num2str(ret1) '.mat']);
    end
    for ret2=1:2
        behav.ret2{id,ret2}  = load([paths.behav IDs{id} '_ret2' num2str(ret2) '.mat']);
    end
    
end


%% analyse

% basic accuracy

% encoding accuracy
encdat=[];
for id=1:length(IDs)
    
    % gather data in one variable
    encdat{id,1}.accuracy = [behav.enc{id,1}.dat.encoding.results.accuracy;behav.enc{id,2}.dat.encoding.results.accuracy];
    encdat{id,1}.rt       = [behav.enc{id,1}.dat.encoding.results.rt_resp;behav.enc{id,2}.dat.encoding.results.rt_resp];
    
    % compute stats
    encAccuracy(id,1)     = nanmean(encdat{id,1}.accuracy);
    encRT(id,1)           = nanmean(encdat{id,1}.rt);
    
end

% retrieval 1 accuracy
ret1dat=[];
for id=1:length(IDs)
    
    % dather data in one variable
    ret1dat{id,1}.accuracy = [behav.ret1{id,1}.dat.retrieval1.results.accuracy;...
        behav.ret1{id,2}.dat.retrieval1.results.accuracy;...
        behav.ret1{id,3}.dat.retrieval1.results.accuracy;...
        behav.ret1{id,4}.dat.retrieval1.results.accuracy]; % accuracy
    ret1dat{id,1}.rt_obj   = [behav.ret1{id,1}.dat.retrieval1.results.rt_resp;...
        behav.ret1{id,2}.dat.retrieval1.results.rt_resp;...
        behav.ret1{id,3}.dat.retrieval1.results.rt_resp;...
        behav.ret1{id,4}.dat.retrieval1.results.rt_resp]; % rt at the object screen
    ret1dat{id,1}.rt_conf  = [behav.ret1{id,1}.dat.retrieval1.results.rt_conf;...
        behav.ret1{id,2}.dat.retrieval1.results.rt_conf;...
        behav.ret1{id,3}.dat.retrieval1.results.rt_conf;...
        behav.ret1{id,4}.dat.retrieval1.results.rt_conf]; % rt at the confidence
    % categories-> 1=cue1, 2=cue2, 3=control, 0=null
    ret1dat{id,1}.trial    = [behav.ret1{id,1}.dat.retrieval1.config.stim.stimvec;...
        behav.ret1{id,2}.dat.retrieval1.config.stim.stimvec;...
        behav.ret1{id,3}.dat.retrieval1.config.stim.stimvec;...
        behav.ret1{id,4}.dat.retrieval1.config.stim.stimvec]; % trial information
    
    % compute accuracy
    clear tmpacc
    rmind = (ret1dat{id,1}.trial==0 | ret1dat{id,1}.trial==3);
    tmpacc = ret1dat{id,1}.accuracy; tmpacc(rmind)=[];
    ret1AccuracyOverall(id,1) = nanmean(tmpacc==1);
    
    clear tmpacc tmptrl
    tmpacc = ret1dat{id,1}.accuracy; tmpacc(rmind)=[];
    tmpacc = tmpacc==1;
    tmptrl = ret1dat{id,1}.trial; tmptrl(rmind)=[];
    ret1AccuracyCue1(id,1) = mean(tmpacc(tmptrl==1));
    ret1AccuracyCue2(id,1) = mean(tmpacc(tmptrl==2));
    
    
    % compute reaction time
    clear tmprt
    tmprt = ret1dat{id,1}.rt_obj; tmprt(rmind)=[];
    nanind = isnan(tmprt);
    tmprt(nanind)=[];
    ret1RTOveralll(id,1) = mean(tmprt);
    tmptrl(nanind)=[];
    ret1RTCue1(id,1) = mean(tmprt(tmptrl==1));
    ret1RTCue2(id,1) = mean(tmprt(tmptrl==2));

end


% retrieval 2 accuracy
ret2dat=[];
for id=1:length(IDs)
    
    % dather data in one variable
    ret2dat{id,1}.accuracy = [behav.ret2{id,1}.dat.retrieval2.results.accuracy;...
        behav.ret2{id,2}.dat.retrieval2.results.accuracy]; % accuracy
    ret2dat{id,1}.rt_obj   = [behav.ret2{id,1}.dat.retrieval2.results.rt_resp;...
        behav.ret2{id,2}.dat.retrieval2.results.rt_resp]; % rt at the object screen
    ret2dat{id,1}.rt_conf  = [behav.ret2{id,1}.dat.retrieval2.results.rt_conf;...
        behav.ret2{id,2}.dat.retrieval2.results.rt_conf]; % rt at the confidence
    % categories-> 1=cue1, 2=cue2, 3=control, 0=null
    ret2dat{id,1}.trial    = [behav.ret2{id,1}.dat.retrieval2.config.stim.stimvec;...
        behav.ret2{id,2}.dat.retrieval2.config.stim.stimvec]; % trial information
    
    % compute accuracy
    clear tmpacc
    rmind = (ret2dat{id,1}.trial==0 | ret2dat{id,1}.trial==3);
    tmpacc = ret2dat{id,1}.accuracy; tmpacc(rmind)=[];
    ret2AccuracyOverall(id,1) = nanmean(tmpacc==1);
    
    clear tmpacc tmptrl
    tmpacc = ret2dat{id,1}.accuracy; tmpacc(rmind)=[];
    tmpacc = tmpacc==1;
    tmptrl = ret2dat{id,1}.trial; tmptrl(rmind)=[];
    ret2AccuracyCue1(id,1) = mean(tmpacc(tmptrl==1));
    ret2AccuracyCue2(id,1) = mean(tmpacc(tmptrl==2));
    
    
    % compute reaction time
    clear tmprt
    tmprt = ret2dat{id,1}.rt_obj; tmprt(rmind)=[];
    nanind = isnan(tmprt);
    tmprt(nanind)=[];
    ret2RTOveralll(id,1) = mean(tmprt);
    tmptrl(nanind)=[];
    ret2RTCue1(id,1) = mean(tmprt(tmptrl==1));
    ret2RTCue2(id,1) = mean(tmprt(tmptrl==2));

end


%% visualise

close all

% retrieval 1 accuracy
figure;
subplot(2,1,1);
clear X Y
X = categorical({'1st Cue','2nd Cue'});
X = reordercats(X,{'1st Cue','2nd Cue'});
Y = [ret1AccuracyCue1 ret1AccuracyCue2];
barlab=[ret1AccuracyCue1; ret1AccuracyCue2];
bar(X,Y)
% text([1:length(barlab)], barlab', num2str(barlab','%0.1f'),'HorizontalAlignment','center','VerticalAlignment','bottom')
% text(1:length(barlab),barlab,num2str(barlab'),'vert','bottom','horiz','center'); 
title('retrieval 1 accuracy, single subjects','FontSize',30)

subplot(2,1,2);
clear X Y
X = categorical({'1st Cue','2nd Cue'});
X = reordercats(X,{'1st Cue','2nd Cue'});
Y = [mean(ret1AccuracyCue1) mean(ret1AccuracyCue2)];
bar(X,Y)
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center'); 
ylim([0 1])
title('retrieval 1 accuracy, mean','FontSize',30)



% retrieval 2 accuracy
figure;
subplot(2,1,1);
clear X Y
X = categorical({'1st Cue','2nd Cue'});
X = reordercats(X,{'1st Cue','2nd Cue'});
Y = [ret2AccuracyCue1 ret2AccuracyCue2];
bar(X,Y)
title('retrieval 2 accuracy, single subjects','FontSize',30)

subplot(2,1,2);
clear X Y
X = categorical({'1st Cue','2nd Cue'});
X = reordercats(X,{'1st Cue','2nd Cue'});
Y = [mean(ret2AccuracyCue1) mean(ret2AccuracyCue2)];
bar(X,Y)
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center'); 
ylim([0 1])
title('retrieval 2 accuracy, mean','FontSize',30)

