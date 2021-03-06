%% Load data

%% Compute likelihood ratios for firing order words
P = M.P;
W = W_Post;
words = words_Post;
clear M W_Post words_Post

%% Remove words less than 2 letters
words = words(cellfun('length',words)>1);
P0 = burst_independent_model(W);

L = zeros(size(words));
for iword = 1:length(words)
    t_prob = sequence_prob(P, words{iword});
    n_prob = sequence_prob(P0, words{iword});
    L(iword) = t_prob / n_prob;
end

likelihood_ratios = [];
likelihood_ratios.fo_words = words;
likelihood_ratios.L = L;

likelihoodProportion = sum(L>1)/length(L);
%%
% cell2mat(cellfun(@(x) [sum(x.L > 1), length(x.L)], likelihood_ratios,...
%     'UniformOutput', false))



