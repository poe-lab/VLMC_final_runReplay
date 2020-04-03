%% Compute likelihood ratios for each CV round

nfold = 10;
d = 2;
alg = 'PST';

RUN_gap_subsequence_cv_data = cell(4, 1);
for rat = 1:4
    rat,
    
    % Compute probabilities
    W = sleep_W{rat};
    
    % Compute folds
    L = length(W);
    fold_size = floor(L / nfold);
    fold = cell(nfold, 1);
    for ifold = 1:nfold
        fold{ifold} = (ifold - 1) * fold_size + (1:fold_size);
    end
    
    % Fit model on each fold
    cv_P = cell(nfold, 1);
    cv_P0 = cell(nfold, 1);
    for ifold = 1:nfold
        train_W = [W(cell2mat(fold(1:(ifold-1))')), zeros(1, d),...
            W(cell2mat(fold((ifold+1):end)'))];
        
        M = word2markovchain(train_W, nfold, d, alg);
        
        cv_P{ifold} = M.P;
        cv_P0{ifold} = burst_independent_model(train_W);
    end
    
    % Words with gaps
    seq = RUN_seq{rat}';
    fo_words = firing_order_words(seq);
    
    % Get cross validation values for log-likelihood ratio
    fit_prob = zeros(length(fo_words), nfold);
    null_prob = zeros(length(fo_words), nfold);
    log_lik = zeros(length(fo_words), nfold);
    for iword = 1:length(fo_words)
        for ifold = 1:nfold
            fit_prob(iword, ifold) = sequence_prob(cv_P{ifold}, fo_words{iword});
            null_prob(iword, ifold) = sequence_prob(cv_P0{ifold}, fo_words{iword});
        end
        log_lik = log2(fit_prob ./ null_prob);
    end
    
    % Z-score
    mn = mean(log_lik, 2);
    sd = std(log_lik, 0, 2);
    z = mn ./ sd;
    
    % Output data structure
    RUN_gap_subsequence_cv_data{rat} = [];
    RUN_gap_subsequence_cv_data{rat}.fo_words = fo_words;
    RUN_gap_subsequence_cv_data{rat}.fit_prob = fit_prob;
    RUN_gap_subsequence_cv_data{rat}.null_prob = null_prob;
    RUN_gap_subsequence_cv_data{rat}.cv_P = cv_P;
    RUN_gap_subsequence_cv_data{rat}.cv_P0 = cv_P0;
    RUN_gap_subsequence_cv_data{rat}.log_lik = log_lik;
    RUN_gap_subsequence_cv_data{rat}.z = z;
    
end

%% Save data

save('RUN_gap_subsequence_cv_data.mat', 'RUN_gap_subsequence_cv_data')

%% Box plot

rat = 1;

log_lik = RUN_gap_subsequence_cv_data{rat}.log_lik;

med_cv = median(log_lik, 2);

%[~, idx] = sort(med_cv);
idx = 1:length(med_cv);

figure;
hold on;
boxplot(log_lik(idx, :)', 'orientation', 'horizontal')
vline(0, 'k');
boxplot(log_lik(idx, :)', 'orientation', 'horizontal')

1 - normcdf(RUN_gap_subsequence_cv_data{rat}.z),

%% Print stats to file

sig_fo_word_stats = ...
    [{'Rat', 'Num > 0', 'Num Sig', 'Longest Sig?', 'Total words'};...
    cell(4, 5)];
for rat = 1:4
    
    fo_words = RUN_gap_subsequence_cv_data{rat}.fo_words;
    log_lik = RUN_gap_subsequence_cv_data{rat}.log_lik;
    
    mn = mean(log_lik, 2);
    sd = std(log_lik, [], 2);
    
    stand_effect_size = mn ./ sd;
    sig_thresh = norminv(1 - 0.05 / length(stand_effect_size));
    
    % Which rat?
    sig_fo_word_stats{rat + 1, 1} = num2str(rat); 
    
    % Number of words with mean log likelihood > 0
    sig_fo_word_stats{rat + 1, 2} = num2str(sum(mn > 0));
    
    % Number of words with significant enrichment
    sig_fo_word_stats{rat + 1, 3} = num2str(sum(stand_effect_size > sig_thresh));
    
    % Is the longest word significant
    sig_fo_word_stats{rat + 1, 4} = 'no';
    if stand_effect_size(end) > sig_thresh
        sig_fo_word_stats{rat + 1, 4} = 'yes';
    end
    
    % Number of words total
    sig_fo_word_stats{rat + 1, 5} = num2str(length(fo_words));
    
end

sig_fo_stats_file = 'sig_fo_stats.txt';
cell2file(sig_fo_stats_file, sig_fo_word_stats);

%% Plot figures

for rat = 1:4
    
    fo_words = RUN_gap_subsequence_cv_data{rat}.fo_words;
    log_lik = RUN_gap_subsequence_cv_data{rat}.log_lik;
    
    mn = mean(log_lik, 2);
    sd = std(log_lik, [], 2);
    
    stand_effect_size = mn ./ sd;
    sig_thresh = norminv(1 - 0.05 / length(stand_effect_size));

    fig = figure;
    hist(stand_effect_size, 20)
    h = vline(sig_thresh, 'r');
    set(h, 'LineWidth', 2)
    set(gca, 'xlim', [-10, 40], 'FontSize', 25, 'LineWidth', 2)

    fig_title = ['Z-scores_distribution_firing_order_words_w_gaps_rat_',...
        num2str(rat)];

    print(fig, '-dtiff', '-r600', [fig_title, '.tif'])
    close(fig)

end




