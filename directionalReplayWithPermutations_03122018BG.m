function directionalReplayWithPermutations_03122018BG(BH_FDR, EschenkoSpindle)%(session,refs)

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
numPermutations = 500;
numWords = length(BH_FDR.wordState);
pVal = cell(numWords,1);
rho = cell(numWords,1);
syllableIdxAll = cell(numWords,1);
replayDirection = cell(numWords,1);
parfor m = 1:numWords
    numSizeSequences = BH_FDR.wordLength(m) - (minLetters - 1); % # of different lengths of syllables from 5:n (n = word length)
%     subSeq = cell(factorial(numSizeSequences,1)); % # of combinations of syllables
    numSyllables = sum(1:numSizeSequences);
    
    % Define output variables with preallocated space:
    syllableIdx = zeros(numSyllables, 2);
    Rs = zeros(numSyllables,1); % Matrix with the correlations (for first spike, mean spike time, median spike)
    Ps = zeros(numSyllables,1); % Matrix with the p-values
    sigTargIdx = zeros(numSyllables,1); % Matrix with the p-values
    %% Create Word permutation matrix
    permutWord = zeros(BH_FDR.wordLength(m), numPermutations);                         % Allocate memory for the permutations (stacked in columns)
    for n = 1:numPermutations                                                         % For each repetition
        z = randperm(BH_FDR.wordLength(m));                                   % Randomly permute the ENTRIES of the cell sequence
        permutWord(:,n) = BH_FDR.words{m}(z)';                                   % Store the shuffled cells
    end

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
           
            %% Create all forward sequences with circular shifting:
            sortSyllable = sort(syllable);
            numShifts = size(sortSyllable,1) - 1;
            forwardSeq = sortSyllable;
            for l = 1:numShifts
                forwardSeq = [forwardSeq circshift(sortSyllable, l)];
            end

            %% Calculate correlation between syllable & all circ shifted forward sequences:
            [rs ,p] = corr(syllable, forwardSeq, 'type','Spearman'); % Spearman is nonparametric
            [Ps(syllableNum),b] = min(p);
            Rs(syllableNum) = rs(b);
            highCorrSeq = forwardSeq(:,b);            
            
            %% Compute correlations of permutations:
            if Ps(syllableNum) < 0.05
                permutRho = zeros(numPermutations,1);  % Matrix with the correlations
                for o = 1:numPermutations  % For each permutation
                    permutRho(o) = corr(permutWord(startSeq:endSeq,o), highCorrSeq, 'type','Spearman'); % Get the correlation of shuffled and the forward sequence
                end

                %% COMPUTE SIGNIFICANCE AND GET DIRECTION OF REPLAY
                % Checks if the correlation of each sequence is significant (above or below
                % 95% of the random permutations correlations)
                Rless = sum(permutRho <= Rs(syllableNum));  % Find how many permutations of the replay have less correlation with the prototype sequence than the original sequence
                Rless = Rless/numPermutations*100;   % Turn it to percentage
                Rmore = sum(permutRho >= Rs(syllableNum));  % Same for permutations with higher correlation than the original
                Rmore = Rmore/numPermutations*100;
                if Rless >= 95   % If there are more than 95% permuations with correlation less than the original sequence
                    sigTargIdx(syllableNum) = 1;  % Mark this replay as a forward replay
                elseif Rmore >= 95  % If there are more than 95% permutations with correlation higher than the original sequence
                    sigTargIdx(syllableNum) = -1;  % Mark this replay as a reverse replay
                end
            end
        end       
    end
    
    %% Keep syllables with p-value < 0.05
    targIdx = sigTargIdx ~= 0;
    pVal{m, 1} = Ps(targIdx);
    rho{m, 1} = Rs(targIdx);
    syllableIdxAll{m, 1} = syllableIdx(targIdx, :);
    replayDirection{m, 1} = sigTargIdx(targIdx);
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
replay.direction = replayDirection(keepData);

%% Find words in spindles
[replay.inSpindle, spindleWordsReplay] = wordsInSpindles_02262018(replay.wordsTimeBins, EschenkoSpindle);


save(['Replays',num2str(refshank),'_',num2str(refchannel),'.mat'],'Event_times','Event_indexes',...
    'Sequences_f','Sequences_r','Sequences_noreps','periripdur');