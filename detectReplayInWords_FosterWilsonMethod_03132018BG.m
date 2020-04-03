function detectReplayInWords_FosterWilsonMethod_03132018BG(BH_FDR, EschenkoSpindle,PC_post)
% This function analyses replay in the word sequences identified as
% significant in the VLMC analyses.
%--------------------------------------------------------------------------
% Created by Brooks A. Gross
% VLMC analyses adapted from Matt Mahoney
% Correlation analyses adapted from Jiannis Taxiditis 
%--------------------------------------------------------------------------

% repet = 500; % Number of random permutations for each replay
% 
% hbins = -1 : 0.05 : 1; % Bins for correlation histograms
% rN = size(targetBins,1);
% disp(' ')
% disp([num2str(rN),' bins in total'])
% disp(' ')

%% Run analysis for each significant word
minLetters = 5;
numCells =length(PC_post);
minCells = max([minLetters ceil(numCells/3)]);
numWords = length(BH_FDR.wordState);
pVal = cell(numWords,1);
rho = cell(numWords,1);
spikeSeq = cell(numWords,1);
fwdSeqHiCorr = cell(numWords,1);

parfor m = 1:numWords
    %% Get all spike bursts for word:
    startTime = BH_FDR.wordsTimeBins(m,1);
    stopTime =  BH_FDR.wordsTimeBins(m,2);
    max_isi = 0.05;
    burstISI = 0.016;
    spikesInWord = [];
    spikeBursts = [];
    
    for c = 1:numCells % For each place cell
        spikes = PC_post{c}; % Get all its spikes
        spikes = spikes(spikes >= startTime & spikes <= stopTime); % Keep all its spikes that fall within the target interval
        
        % Combine all spikes from all place cells into master spike train:
        if ~isempty(spikes) % If there are any spikes
            spikesInWord = [spikesInWord; [spikes c*ones(size(spikes,1),1)]];
        end
        % Replay coding method: Cleaned sequence by keeping 1st spike time of each burst
        if length(spikes) >= 1
            isiUnit = diff(spikes);
            keep_idx = find([1; isiUnit > burstISI] > 0);
            spikes =spikes(keep_idx);
            spikeBursts = [spikeBursts; [spikes c*ones(size(spikes,1),1)]];
        end
    end
    [~, firingOrder] = sort(spikesInWord(:,1));
    spikesInWord = spikesInWord(firingOrder,:);
    
    [~, firingOrder] = sort(spikeBursts(:,1));
    spikeBursts = spikeBursts(firingOrder,:);
    
    %% FOSTER and WILSON method:
    % Separate into bins where ISI is greater than 50 msec
    lengthSpikesInWord = size(spikesInWord,1);
    isi = diff(spikesInWord(:,1));
    binBreaks = find(isi > max_isi);
    if isempty(binBreaks)
        replayBins = [1 lengthSpikesInWord];
    else    
        replayBins = [[1;binBreaks+1], [binBreaks;lengthSpikesInWord]];
    end
    
    % Remove bins longer than 500 msec
    replayDuration = spikesInWord(replayBins(:,2)) - spikesInWord(replayBins(:,1));
    replayBins = replayBins(replayDuration<.5,:);
    
    % Remove bins where less than max(5 or 1/3 of the place cells) fire at least once
    numUniqueCells = [];
    for i = 1:size(replayBins,1) 
        numUniqueCells(i) = length(unique(spikesInWord(replayBins(i,1):replayBins(i,2) , 2)));
    end
    replayBins = replayBins(numUniqueCells>=minCells,:);
    
    % Get bin time boundaries:
    [numReplayBins, binCol] = size(replayBins);
    replayBinTS = zeros(numReplayBins, binCol);
    for i = 1:numReplayBins
        for j = 1:binCol
            replayBinTS(i,j) = spikesInWord(replayBins(i,j),1);    
        end
    end
    
    % Define output variables with preallocated space:
    Rs = zeros(numReplayBins,1); % Matrix with the correlations (for first spike, mean spike time, median spike)
    Ps = zeros(numReplayBins,1); % Matrix with the p-values
    highCorrSeq = cell(numReplayBins,1); % Highest correlated forward replay sequences
    
    % Define each possible combination of syllables of >= 5 letters and run
    % correlation analyses to determine if significant in forward or
    % reverse direction:
    for i= 1:numReplayBins            
        tempReplay = spikeBursts((spikeBursts(:,1) >= replayBinTS(i,1) &...
            spikeBursts(:,1) <= replayBinTS(i,2)),:);
        
        %% Create all forward sequences with circular shifting
        sortReplaySeq = sort(tempReplay(:,2));
        numShifts = size(sortReplaySeq,1) - 1;
        forwardSeq = sortReplaySeq;
        for j = 1:numShifts
            forwardSeq = [forwardSeq circshift(sortReplaySeq,j)];
        end
        
        %% Calculate correlation between syllable & all circ shifted forward sequences:
        [rs ,p] = corr(tempReplay(:,2), forwardSeq, 'type','Spearman'); % Spearman is nonparametric
        [Ps(i),b] = min(p);
        Rs(i) = rs(b);
        highCorrSeq{i,1} = forwardSeq(:,b);
%         clear b rs p syllable syllableNum forwardSeq  
    end
    
    %% Keep syllables with p-value < 0.05
    targIdx = Ps < 0.05;
    pVal{m,1} = Ps(targIdx);
    rho{m,1} = Rs(targIdx);
    fwdSeqHiCorr{m,1} = highCorrSeq(targIdx);
    spikeSeq{m,1} = spikeBursts;
%     clear targIdx Ps Rs syllableIdx
end

%% Remove all significant words that did not have any directional replay:
cellfun(@(x) ~isempty(x), pVal);
keepData = cellfun(@(x) ~isempty(x), pVal);
replay.pVal = pVal(keepData);
replay.rho = rho(keepData);
replay.fwdSeqHiCorr = fwdSeqHiCorr(keepData);
replay.spikeSeq = spikeSeq(keepData);
replay.words = BH_FDR.words(keepData);
replay.wordsTimeBins = BH_FDR.wordsTimeBins(keepData,:);
replay.wordState = BH_FDR.wordState(keepData);
replay.wordLength = BH_FDR.wordLength(keepData);

%% Find words in spindles
[replay.inSpindle, spindleWordsReplay] = wordsInSpindles_02262018(replay.wordsTimeBins, EschenkoSpindle);
%% SAVE
save(['Replays',num2str(refshank),'_',num2str(refchannel),'.mat'],'Event_times','Event_indexes',...
    'Sequences_f','Sequences_r','Sequences_noreps','periripdur');