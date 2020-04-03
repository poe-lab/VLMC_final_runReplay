function [replay, spindleWordsReplay] = detectReplayInWordsSyll_FosterWilsonMethod_06072018BG(BH_FDR, EschenkoSpindle,PC_post)
% This function analyzes replay in the word sequences identified as
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
minCells = 5; %max([minLetters ceil(numCells/3)]);
numWords = length(BH_FDR.wordState);
pVal = cell(numWords,1);
rho = cell(numWords,1);
spikeSeq = cell(numWords,1);
fwdSeqHiCorr = cell(numWords,1);
syllableIdxAll = cell(numWords,1);
%% Use 'parfor' parallell computing command for the following 'For' loop:
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
    highCorrSeqRB = cell(numReplayBins,1); % Highest correlated forward replay sequences
    pValRB = cell(numReplayBins,1);
    rhoRB = cell(numReplayBins,1);
    spikeSeqRB = cell(numReplayBins,1);
    syllableIdxRB = cell(numReplayBins,1);
    
    % Define each possible combination of syllables of >= 5 letters and run
    % correlation analyses to determine if significant in forward or
    % reverse direction:
    for i= 1:numReplayBins            
        tempReplay = spikeBursts((spikeBursts(:,1) >= replayBinTS(i,1) &...
            spikeBursts(:,1) <= replayBinTS(i,2)),:);
        
        lengthReplay = size(tempReplay,1); 
        numSizeSequences = lengthReplay - (minLetters - 1); % # of different lengths of syllables from 5:n (n = word length)
        numSyllables = sum(1:numSizeSequences);
        
        % Define output variables with preallocated space:
        syllableIdx = zeros(numSyllables, 2);
        Rs = zeros(numSyllables,1); % Matrix with the correlations (for first spike, mean spike time, median spike)
        Ps = -1*ones(numSyllables,1); % Matrix with the p-values
        highCorrSyll = cell(numSyllables,1);
        % Define each possible combination of syllables of >= 5 letters and run
        % correlation analyses to determine if significant in forward or
        % reverse direction:
        for j= 1:numSizeSequences
            sizeSubSeq = lengthReplay - (j - 1);
        
            for k = 1:j
                startSeq = k;
                endSeq = k - 1 + sizeSubSeq;
                syllableNum = sum(1:j-1) + k;
                syllableIdx(syllableNum, :) = [startSeq endSeq];
                syllable = tempReplay(startSeq:endSeq,2);
                if isequal(syllable(1), syllable(2)) || isequal(syllable(end), syllable(end-1))
                    % Skip syllables with 2 or more of the same letter
                    % (place cell) at the beginning or end of the sequence
                else
                    % Condense sequence: Remove adjacent repeats of the
                    % same letter (place cell) from the sequence before
                    % calculating the correlation:
                    condensedSyll = syllable(logical([1; diff(syllable) ~= 0]));
                    if length(condensedSyll) >= minLetters
                        %% Create all forward sequences with circular shifting
                        sortSyllable = sort(condensedSyll);
                        numFwdSeq = size(sortSyllable,1);
                        forwardSeq = [];
                        for l = 1:numFwdSeq
                            circShftSeq = circshift(sortSyllable, l-1);   
                            shftNoRepeat = circShftSeq(logical([1; diff(circShftSeq) ~= 0]));
                            if isequal(circShftSeq, shftNoRepeat)
                                forwardSeq = [forwardSeq circShftSeq];
                            end
                        end
                        
                        
%                         numShifts = size(sortSyllable,1) - 1;
%                         forwardSeq = sortSyllable;
%                         
%                         for l = 1:numShifts
%                             circShftSeq = circshift(sortSyllable, l);   
%                             shftNoRepeat = circShftSeq(logical([1; diff(circShftSeq) ~= 0]));
%                             if isequal(circShftSeq, shftNoRepeat)
%                                 forwardSeq = [forwardSeq circShftSeq];
%                             end
%                         end

                        %% Calculate correlation between syllable & all circ shifted forward sequences:
                        if ~isempty(forwardSeq)
                            [rs ,p] = corr(condensedSyll, forwardSeq, 'type','Spearman'); % Spearman is nonparametric
                            [Ps(syllableNum),b] = min(p);
                            Rs(syllableNum) = rs(b);
                            highCorrSyll{syllableNum,1} = forwardSeq(:,b);
                        end
                    end
                end    
            end
            
        end
        %% Keep syllables with p-value < 0.05
        targIdx = Ps < 0.05 & Ps>=0;
        pValRB{i,1} = Ps(targIdx);
        rhoRB{i,1} = Rs(targIdx);
        syllableIdxRB{i,1} = syllableIdx(targIdx,:);
        highCorrSeqRB{i,1} = highCorrSyll(targIdx);
        spikeSeqRB{i,1} = tempReplay;
    end
    keepData = cellfun(@(x) ~isempty(x), pValRB);
    pVal{m,1} = pValRB(keepData);
    rho{m,1} = rhoRB(keepData);
    syllableIdxAll{m,1} = syllableIdxRB(keepData);
    fwdSeqHiCorr{m,1} = highCorrSeqRB(keepData);
    spikeSeq{m,1} = spikeSeqRB(keepData); 
end

%% Remove all significant words that did not have any directional replay:
keepData = cellfun(@(x) ~isempty(x), pVal);
replay.pVal = pVal(keepData);
replay.rho = rho(keepData);
replay.syllableIdxAll = syllableIdxAll(keepData);
replay.fwdSeqHiCorr = fwdSeqHiCorr(keepData);
replay.spikeSeq = spikeSeq(keepData);
replay.words = BH_FDR.words(keepData);
replay.wordsTimeBins = BH_FDR.wordsTimeBins(keepData,:);
replay.wordState = BH_FDR.wordState(keepData);
replay.wordLength = BH_FDR.wordLength(keepData);

%% Find words in spindles
[replay.inSpindle, spindleWordsReplay] = wordsInSpindles_02262018(replay.wordsTimeBins, EschenkoSpindle);
% %% SAVE
% save(['Replays',num2str(refshank),'_',num2str(refchannel),'.mat'],'Event_times','Event_indexes',...
%     'Sequences_f','Sequences_r','Sequences_noreps','periripdur');