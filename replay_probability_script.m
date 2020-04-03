%% Load data

cd ~/Desktop/VLMC_final/Figures/RunReplay/

load ../RunRasters/RUN_sequences.mat
load ../Words/Sleep_words_and_sequences.mat
load ../VlmcFitting/Second_order_Markov_chains.mat

addpath('~/Desktop/VLMC_final');

%% RUN sequence probabilities

RUN_seq_prob = zeros(4, 1);
for rat = 1:4
    
    seq = RUN_seq{rat}';
    P = second_order_models{rat};
    RUN_seq_prob(rat) = sequence_prob(P, seq);
    
end

%% Null distribution for probabilities

nsamp = 20;

nfold = 10;
d = 0:5;
alg = 'PST';

RUN_seq_prob_null = zeros(4, nsamp);
for rat = 1:4
    
    seq = RUN_seq{rat}';
    P = second_order_models{rat};
    W = sleep_W{rat};
    
    for isamp = 1:nsamp
        
        [rat, isamp],
        
        null_P = scramble_sequence_markov_null(W);
        
        RUN_seq_prob_null(rat, isamp) = sequence_prob(null_P, seq);
        
    end
end

%% Save

save('Null_distributions_for_RUN_sequence_probability.mat', 'RUN_seq_prob_null')

%%

emp_pval = zeros(4, 1);
the_pval = zeros(4, 1);
for rat = 1:4
    
    curr_val = RUN_seq_prob(rat);
    curr_null = RUN_seq_prob_null(rat, :);
    
    emp_pval(rat) = sum(curr_null > curr_val) / nsamp;
    
    mn = mean(curr_null);
    sd = std(curr_null);
    
    the_pval(rat) = 1 - normcdf((curr_val - mn) / sd);
    
end







