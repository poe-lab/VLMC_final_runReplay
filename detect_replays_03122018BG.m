function detect_replays_03122018BG(BH_FDR, EschenkoSpindle)%(session,refs)

%--------------------------------------------------------------------------
% Created by Jiannis Taxidis, Caltech, USA, February 2013
% Adapted for VLMC significant words
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
numWords = length(BH_FDR.wordState);
pVal = cell(numWords,1);
rho = cell(numWords,1);
syllableIdxAll = cell(numWords,1);
parfor m = 1:numWords
    numSizeSequences = BH_FDR.wordLength(m) - (minLetters - 1); % # of different lengths of syllables from 5:n (n = word length)
%     subSeq = cell(factorial(numSizeSequences,1)); % # of combinations of syllables
    numSyllables = sum(1:numSizeSequences);
    % Define output variables with preallocated space:
    syllableIdx = zeros(numSyllables, 2);
    Rs = zeros(numSyllables,1); % Matrix with the correlations (for first spike, mean spike time, median spike)
    Ps = zeros(numSyllables,1); % Matrix with the p-values
    
    % Define each possible combination of syllables of >= 5 letters and run
    % correlation analyses to determine if significant in forward or
    % reverse direction:
    for j= 1:numSizeSequences
        sizeSubSeq = BH_FDR.wordLength(m) - (j - 1);
               
        for k = 1:j
            startSeq = k;
            endSeq = k - 1 + sizeSubSeq;
            syllableNum = sum(1:j-1) + k;
            syllableIdx(syllableNum, :) = [startSeq endSeq];
            syllable = BH_FDR.words{m}(startSeq:endSeq)';
%             clear startSeq endSeq
            
            %% Create all forward sequences with circular shifting
            sortSyllable = sort(syllable);
            numShifts = size(sortSyllable,1) - 1;
            forwardSeq = sortSyllable;
            for l = 1:numShifts
                forwardSeq = [forwardSeq circshift(sortSyllable, l)];
            end
%             clear sortSyllable numShifts

%             highCorrSeq = cell(numSeq,1);
            %% Calculate correlation between syllable & all circ shifted forward sequences:
            [rs ,p] = corr(syllable, forwardSeq, 'type','Spearman'); % Spearman is nonparametric
            [Ps(syllableNum),b] = min(p);
            Rs(syllableNum) = rs(b);
%             highCorrSeq{r} = forwardSeq{r}(:,b);
%             clear b rs p syllable syllableNum forwardSeq
        end       
    end
    
    %% Keep syllables with p-value < 0.05
    targIdx = Ps < 0.05;
    pVal{m, 1} = Ps(targIdx);
    rho{m, 1} = Rs(targIdx);
    syllableIdxAll{m, 1} = syllableIdx(targIdx, :);
%     clear targIdx Ps Rs syllableIdx
end

%% Remove all significant words that did not have any directional replay:
cellfun(@(x) ~isempty(x), pVal);
keepData = cellfun(@(x) ~isempty(x), pVal);
replay.pVal = pVal(keepData);
replay.rho = rho(keepData);
replay.syllableIdxAll = syllableIdxAll(keepData);
replay.words = BH_FDR.words(keepData);
replay.wordsTimeBins = BH_FDR.wordsTimeBins(keepData,:);
replay.wordState = BH_FDR.wordState(keepData);
replay.wordLength = BH_FDR.wordLength(keepData);

%% Find words in spindles
[replay.inSpindle, spindleWordsReplay] = wordsInSpindles_02262018(replay.wordsTimeBins, EschenkoSpindle);
%% SAVE
save(['Replays',num2str(refshank),'_',num2str(refchannel),'.mat'],'Event_times','Event_indexes',...
    'Sequences_f','Sequences_r','Sequences_noreps','periripdur');