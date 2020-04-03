%% Load data

load RUN_gap_subsequence_cv_data.mat
load ../VlmcFitting/Second_order_Markov_chains.mat
load ../Words/Sleep_words_and_sequences.mat

%% Compute likelihood ratios for firing order words

RUN_order_likelihood_ratios = cell(4, 1);
for rat = 1:4
    rat,
    
    P = second_order_models{rat};
    
    W = sleep_W{rat};
    P0 = burst_independent_model(W);
    
    fo_words = RUN_gap_subsequence_cv_data{rat}.fo_words;
    
    L = zeros(size(fo_words));
    for iword = 1:length(fo_words)
        
        t_prob = sequence_prob(P, fo_words{iword});
        n_prob = sequence_prob(P0, fo_words{iword});
        
        L(iword) = t_prob / n_prob;
        
    end
    
    RUN_order_likelihood_ratios{rat} = [];
    RUN_order_likelihood_ratios{rat}.fo_words = fo_words;
    RUN_order_likelihood_ratios{rat}.L = L;
    
end

%%

cell2mat(cellfun(@(x) [sum(x.L > 1), length(x.L)], RUN_order_likelihood_ratios,...
    'UniformOutput', false))



