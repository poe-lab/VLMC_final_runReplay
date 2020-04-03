function SignificantSequences_REM_03192018(words, wordState, wordsTimeBins, M, null_P)

% Get Markov chains
P = M.P;
clear M

%% Select words of target state(s):
targetLogic = wordState == 3;
words = words(targetLogic);
wordsTimeBins = wordsTimeBins(targetLogic,:);
wordState = wordState(targetLogic);
clear targetLogic

%% Remove words less than _ letters
minLetters = 5;
targetLogic = cellfun('length',words)>= minLetters;
words = words(targetLogic);
wordsTimeBins = wordsTimeBins(targetLogic,:);
wordState = wordState(targetLogic);
clear targetLogic

%% Compute word probabilities
word_prob = sequence_prob_BGmod(P, words);

%% Null probability distributions
null_dist = [];
parfor isamp = 1:length(null_P)    
    % Compute probabilities
    word_null_prob = sequence_prob_BGmod(null_P{isamp}, words);
    null_dist = [null_dist, word_null_prob];
end

%% Z-score
true_z = zeros(size(word_prob));
for i = 1:length(word_prob)
    mn = mean(null_dist(i, :));
    sd = std(null_dist(i, :));
    true_z(i) = (word_prob(i) - mn) / sd;
end
pval = 1 - normcdf(true_z);

%% Compute False Discovery Rate - Benjamani-Hochberg Procedure
% The Benjamini–Hochberg procedure (BH step-up procedure) controls the FDR
% at level alpha. It works as follows:
% 1) For a given alpha, find the largest k such that P(k) <= (k/m)*alpha.
% 2) Reject the null hypothesis (i.e., declare discoveries) for all H(i) for i = 1,...,k.
alpha = 0.05;
numTests = length(pval); % m
pThreshold = ((1:1:numTests)/numTests * alpha);
[sort_pVal, idx_pVal]  = sort(pval);
targLogic = sort_pVal <= pThreshold';
targEnd = find(targLogic, 1, 'last');

BH_FDR.words = words(idx_pVal(1:targEnd));
BH_FDR.wordsTimeBins = wordsTimeBins(idx_pVal(1:targEnd),:);
BH_FDR.wordState = wordState(idx_pVal(1:targEnd));

% Sort by time
[~, idx_time]  = sort(BH_FDR.wordsTimeBins(:,1));
BH_FDR.words = BH_FDR.words(idx_time);
BH_FDR.wordsTimeBins = BH_FDR.wordsTimeBins(idx_time,:);
BH_FDR.wordState = BH_FDR.wordState(idx_time);
BH_FDR.wordLength = cellfun('length',BH_FDR.words);







