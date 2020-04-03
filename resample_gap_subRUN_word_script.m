%%

if 0

RUN_gap_subsequence_resampling_data = cell(4, 1);

nsamp = 100;

for rat = 1:4
    % Get Markov chains
    P = second_order_models{rat};
    [s2, s1, P1] = so_markov_moments(P);
    
    F = markov2fut_dist(P);
    
    W = sleep_W{rat};
    alph = unique(W);
    digram = list_kgrams(alph, 2);
    
    seq = RUN_seq{rat}';
    L = length(seq);
    
    % Words with gaps
    fo_words = cell(L, 1);
    for i = 1:(2^L - 1)
        curr_subset = dec2bin(i, L) == '1';
    
        word_size = sum(curr_subset);
    
        fo_words{word_size} = [fo_words{word_size}; seq(curr_subset)];
    end
    fo_words = fo_words(2:end);
    
    % Compute probabilities
    fo_word_prob = [];
    for i = 1:length(fo_words)
        for j = 1:size(fo_words{i}, 1)
            prob = sequence_prob(P, fo_words{i}(j, :));
            fo_word_prob = [fo_word_prob; prob];
        end
    end
    
    % Null distributions
    idx = find(all(digram ~= 0, 2));
    
    null_dist = [];
    for isamp = 1:nsamp
        
        [rat, isamp],
        null_P = scramble_sequence_markov_null(W);
        
        % Compute probabilities
        fo_word_null = [];
        for i = 1:length(fo_words)
            for j = 1:size(fo_words{i}, 1)
                prob = sequence_prob(null_P, fo_words{i}(j, :));
                fo_word_null = [fo_word_null; prob];
            end
        end
        
        null_dist = [null_dist, fo_word_null];
    end
    
    RUN_gap_subsequence_resampling_data{rat} = [];
    RUN_gap_subsequence_resampling_data{rat}.fo_words = fo_words;
    RUN_gap_subsequence_resampling_data{rat}.fo_word_prob = fo_word_prob;
    RUN_gap_subsequence_resampling_data{rat}.null_dist = null_dist;
    
end

%% Z-score

for rat = 1:4
    
    null_dist = RUN_gap_subsequence_resampling_data{rat}.null_dist;
    fo_word_prob = RUN_gap_subsequence_resampling_data{rat}.fo_word_prob;
    
    true_z = zeros(size(fo_word_prob));
    for i = 1:length(fo_word_prob)
        mn = mean(null_dist(i, :));
        sd = std(null_dist(i, :));
        true_z(i) = (fo_word_prob(i) - mn) / sd;
    end
    
    null_z = zscore(null_dist, [], 2);
    
    pval = 1 - normcdf(RUN_gap_subsequence_resampling_data{rat}.true_z);
    
    RUN_gap_subsequence_resampling_data{rat}.pval = pval;
    RUN_gap_subsequence_resampling_data{rat}.null_z = null_z;
    RUN_gap_subsequence_resampling_data{rat}.true_z = true_z;
end

save('RUN_gap_subsequence_resampling_data.mat', 'RUN_gap_subsequence_resampling_data')

end

%% Compute FDR and flatten firing order words cell array

load('RUN_gap_subsequence_resampling_data.mat')

for rat = 1:4
    
    seq = RUN_seq{rat}';
    
    fo_words = RUN_gap_subsequence_resampling_data{rat}.fo_words;
    null_dist = RUN_gap_subsequence_resampling_data{rat}.null_dist;
    fo_word_prob = RUN_gap_subsequence_resampling_data{rat}.fo_word_prob;
    null_z = RUN_gap_subsequence_resampling_data{rat}.null_z;
    true_z = RUN_gap_subsequence_resampling_data{rat}.true_z;
    
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
    
    fo_words = cellfun(...
        @(x) mat2cell(x, ones(size(x, 1), 1), size(x, 2)), fo_words, 'UniformOutput', 0);
    temp_fo_words = [];
    for i = 1:length(fo_words)
        temp_fo_words = [temp_fo_words; fo_words{i}];
    end
    fo_words = temp_fo_words;
    
    fdr = ave_fd ./ (1:length(ave_fd))';
    
    %[sort_z, ave_fd, fdr],
    
    fo_words{idx(1:sum(fdr < 0.1))}
    
    RUN_gap_subsequence_resampling_data{rat}.fdr = fdr;
    RUN_gap_subsequence_resampling_data{rat}.fo_words = fo_words;
    
end

save('RUN_gap_subsequence_resampling_data.mat', 'RUN_gap_subsequence_resampling_data')

