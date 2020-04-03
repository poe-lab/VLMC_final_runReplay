
% Get Markov chains
P = M.P;
clear M

W = W_Post;
words = words_Post;
clear W_Post words_Post

%% Select words of target state(s):
targetLogic = wordState == 2 | wordState == 6;
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

% alph = unique(W);
% digram = list_kgrams(alph, 2);

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

%% Compute FDR - Matt's version
null_z = zscore(null_dist, [], 2);    
% figure;
% hist(mean(null_z))
% vline(mean(true_z), 'r');
%
% disp('Mean p-value')
% sum(mean(null_z) > mean(true_z)) / nsamp,
%
% figure;
% for i = 1:size(null_dist, 1)
%
%     subplot(5, 5, i)
%     hist(null_dist(i, :));
%     vline(fo_word_prob(i), 'r');
%
% end

%     figure;
%     hist(null_z(:), 20)
%     vline(true_z, 'r');

[sort_z, idx]  = sort(true_z, 'descend');

ave_fd = zeros(size(sort_z));
for i = 1:length(sort_z)
    ave_fd(i) = mean(sum(null_z >= sort_z(i)));
end

falseDiscovery.fdr = ave_fd ./ (1:length(ave_fd))';

%[sort_z, ave_fd, fdr],
targEnd = sum(falseDiscovery.fdr < 0.1);
falseDiscovery.words = words(idx(1:targEnd));
falseDiscovery.wordsTimeBins = wordsTimeBins(idx(1:targEnd),:);
falseDiscovery.wordState = wordState(idx(1:targEnd));

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







