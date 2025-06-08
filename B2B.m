%% MIN535 FINAL PROJECT - BRAIN-TO-BRAIN INTERACTION ANALYSIS
% Critical analysis of Tressoldi et al. (2014)
% Final version: 08.06.2025 (Epoch duration set to 4s, robust mean_epoch initialization)
% No Time-Frequency Plots, added Segment-Based Summary Table
clear; clc; close all;

%% =================== SECTION 1: SETTINGS AND CONFIGURATION ===================
Config.fs = 128; % Sampling frequency (Hz)
Config.base_path = 'C:\Users\SmyTp\OneDrive\Masaüstü\Project\1030617\BBI_RawData\BBI_RawData'; % UPDATE YOUR DATA PATH

% Analysis Parameters
Config.fft.nfft = 512;
Config.fft.win = hamming(256);
Config.fft.noverlap = 128;
Config.bands = struct('Delta',[1 4],'Theta',[4 8],'Alpha',[8 13],...
                       'Beta',[13 30],'Gamma',[30 45]);
Config.bnames = fieldnames(Config.bands);

% Segmentation Parameter (as per paper)
Config.segment_len_sec = 4; % Segment length in seconds (4s as per paper)
Config.segment_len_samples = Config.segment_len_sec * Config.fs; % Segment length in samples

% Defining Pairs (Example - Adjust according to your data)
pairs_to_analyze = struct(...
    'sender_id',   {1, 3, 5, 1, 3}, ...
    'receiver_id', {1, 3, 5, 2, 4}, ...
    'type',        {'Paired', 'Paired', 'Paired', 'Unpaired', 'Unpaired'} ...
);

%% =============== SECTION 2: MAIN ANALYSIS LOOP =============
% Expanded results struct to store data for new analyses
results = struct('pair_type', {}, 'psd_corr', {}, 'coh', {}, ...
                 'psd_s_spectrum', {}, 'psd_r_spectrum', {}, 'coh_spectrum', {}, ...
                 'freq_psd', {}, 'freq_coh', {}, ...
                 'mean_s_epoch', {}, 'mean_r_epoch', {}, 'epoch_time', {}); % Removed TF related fields

fprintf('Analysis starting...\n');
for p = 1:length(pairs_to_analyze)
    s_id = pairs_to_analyze(p).sender_id;
    r_id = pairs_to_analyze(p).receiver_id;
    pair_type = pairs_to_analyze(p).type;

    fprintf('Processing: S%d-R%d (%s)\n', s_id, r_id, pair_type);

    % Initialize mean_s_epoch and related variables for each iteration
    % Moved outside try-catch to ensure they are always initialized in scope
    mean_s_epoch = []; 
    mean_r_epoch = []; 
    epoch_time = [];

    try
        % Load Data
        s_path = fullfile(Config.base_path, sprintf('Pair%d', s_id), sprintf('t%d.xlsx', s_id));
        r_path = fullfile(Config.base_path, sprintf('Pair%d', r_id), sprintf('r%d.xlsx', r_id));

        Ts = readtable(s_path, 'VariableNamingRule', 'preserve');
        Tr = readtable(r_path, 'VariableNamingRule', 'preserve');

        eeg_s_raw = table2array(Ts(:, 5:end));
        eeg_r_raw = table2array(Tr(:, 5:end));
        labels_s = Ts{:, 3}; labels_s(labels_s == 2000) = 1;
        labels_r = Tr{:, 3}; labels_r(labels_r == 2000) = 1;

        % Data Dimension Check
        min_len = min(length(labels_s), length(labels_r));
        eeg_s_raw = eeg_s_raw(1:min_len, :);
        eeg_r_raw = eeg_r_raw(1:min_len, :);
        labels = labels_s(1:min_len); % Assume same stimulus labels for both sender and receiver

        % Preprocessing
        eeg_s = preprocess_eeg(eeg_s_raw, Config.fs);
        eeg_r = preprocess_eeg(eeg_r_raw, Config.fs);

        % Apply 4-second segmentation *before* main analyses
        % Only consider data where stimulus is present (labels == 1)
        eeg_s_stim = eeg_s(labels == 1, :);
        eeg_r_stim = eeg_r(labels == 1, :);
        
        [segmented_s_eeg, segmented_r_eeg, num_segments] = segment_eeg_data(eeg_s_stim, eeg_r_stim, Config.segment_len_samples);
        
        if num_segments == 0
            warning('No valid segments found for S%d-R%d. Skipping analysis.', s_id, r_id);
            continue; % Skip to next pair if no valid segments
        end

        % Initialize accumulators for segment-averaged results
        avg_psd_corr_struct = struct();
        avg_coh_struct = struct();
        avg_log_psd_s_all = zeros(Config.fft.nfft/2 + 1, size(eeg_s, 2));
        avg_log_psd_r_all = zeros(Config.fft.nfft/2 + 1, size(eeg_r, 2));
        avg_Cxy_all = zeros(Config.fft.nfft/2 + 1, 1); % Average across channels

        for seg_idx = 1:num_segments
            current_s_segment = segmented_s_eeg(:, :, seg_idx);
            current_r_segment = segmented_r_eeg(:, :, seg_idx);
            
            % Each segment is treated as a 'stimulus' period for spectral analysis
            % labels are not relevant within a segment, so we pass a dummy 'labels_segment'
            labels_segment = ones(size(current_s_segment, 1), 1); 

            % Analyze each segment
            [psd_corr_seg, log_psd_s_seg, log_psd_r_seg, f_psd_seg] = ...
                analyze_psd_correlation(current_s_segment, current_r_segment, labels_segment, Config);
            [coh_seg, Cxy_all_seg, f_coh_seg] = ...
                analyze_coherence(current_s_segment, current_r_segment, labels_segment, Config);
            
            % Accumulate results
            for b_idx = 1:length(Config.bnames)
                band_name = Config.bnames{b_idx};
                if seg_idx == 1
                    avg_psd_corr_struct.(band_name) = psd_corr_seg.(band_name);
                    avg_coh_struct.(band_name) = coh_seg.(band_name);
                else
                    avg_psd_corr_struct.(band_name) = avg_psd_corr_struct.(band_name) + psd_corr_seg.(band_name);
                    avg_coh_struct.(band_name) = avg_coh_struct.(band_name) + coh_seg.(band_name);
                end
            end
            avg_log_psd_s_all = avg_log_psd_s_all + log_psd_s_seg;
            avg_log_psd_r_all = avg_log_psd_r_all + log_psd_r_seg;
            avg_Cxy_all = avg_Cxy_all + Cxy_all_seg;
        end

        % Average accumulated results over segments
        psd_corr_final = avg_psd_corr_struct;
        coh_final = avg_coh_struct;
        for b_idx = 1:length(Config.bnames)
            band_name = Config.bnames{b_idx};
            psd_corr_final.(band_name) = psd_corr_final.(band_name) / num_segments;
            coh_final.(band_name) = coh_final.(band_name) / num_segments;
        end
        log_psd_s_all_final = avg_log_psd_s_all / num_segments;
        log_psd_r_all_final = avg_log_psd_r_all / num_segments;
        Cxy_all_final = avg_Cxy_all / num_segments;
        f_psd = f_psd_seg; % Freq axes are the same for all segments
        f_coh = f_coh_seg;
        
        % Signal Level Data Collection (Average Epochs - now 4 seconds)
        stim_indices = find(labels == 1);
        epoch_length_samples = Config.fs * 4; % <-- Changed epoch duration to 4 seconds
        
        if ~isempty(stim_indices)
            eeg_s_epochs = zeros(length(stim_indices), epoch_length_samples);
            eeg_r_epochs = zeros(length(stim_indices), epoch_length_samples);

            valid_epoch_count = 0;
            for i_stim = 1:length(stim_indices)
                start_idx = stim_indices(i_stim);
                end_idx = min(start_idx + epoch_length_samples - 1, size(eeg_s, 1));
                if (end_idx - start_idx + 1) == epoch_length_samples
                    valid_epoch_count = valid_epoch_count + 1;
                    eeg_s_epochs(valid_epoch_count, :) = mean(eeg_s(start_idx:end_idx, :), 2);
                    eeg_r_epochs(valid_epoch_count, :) = mean(eeg_r(start_idx:end_idx, :), 2);
                end
            end
            
            eeg_s_epochs = eeg_s_epochs(1:valid_epoch_count, :);
            eeg_r_epochs = eeg_r_epochs(1:valid_epoch_count, :);

            mean_s_epoch = mean(eeg_s_epochs, 1);
            mean_r_epoch = mean(eeg_r_epochs, 1);
            epoch_time = (0:epoch_length_samples-1)/Config.fs;
        end

        % Store Results
        results(end+1) = struct(...
            'pair_type', pair_type, ...
            'psd_corr', psd_corr_final, ... % Use segment-averaged results
            'coh', coh_final, ... % Use segment-averaged results
            'psd_s_spectrum', log_psd_s_all_final, ...
            'psd_r_spectrum', log_psd_r_all_final, ...
            'coh_spectrum', Cxy_all_final, ...
            'freq_psd', f_psd, ...
            'freq_coh', f_coh, ...
            'mean_s_epoch', mean_s_epoch, ...
            'mean_r_epoch', mean_r_epoch, ...
            'epoch_time', epoch_time); % Removed TF related fields from results struct

    catch ME
        fprintf('ERROR: S%d-R%d could not be processed: %s\n', s_id, r_id, ME.message);
    end
end


%% ============ SECTION 3: RESULTS ANALYSIS (Visualizations) ============
% Separate Groups
is_paired = strcmp({results.pair_type}, 'Paired');
paired_results = results(is_paired);
unpaired_results = results(~is_paired);

% Metrics
metrics.psd = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
metrics.coh = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};

% Visualization - Boxplots
figure('Name', 'Group Comparison - Band Averages', 'Position', [100, 100, 1200, 500]);
% PSD Correlation Plot
subplot(1,2,1);
hold on;
all_p_data_psd = [];
all_up_data_psd = [];
group_labels_psd = {};
for b = 1:length(metrics.psd)
    band = metrics.psd{b};
    p_data = arrayfun(@(x) x.psd_corr.(band), paired_results);
    up_data = arrayfun(@(x) x.psd_corr.(band), unpaired_results);

    all_p_data_psd = [all_p_data_psd, p_data];
    all_up_data_psd = [all_up_data_psd, up_data];
    group_labels_psd = [group_labels_psd, repmat({[band ' Paired']}, 1, length(p_data)), ...
                          repmat({[band ' Unpaired']}, 1, length(up_data))];

    [~, p] = ttest2(p_data, up_data);
    text(b*2-0.5, max([p_data, up_data])+0.05, sprintf('p=%.3f', p), ...
        'HorizontalAlignment', 'center');
end
boxplot([all_p_data_psd, all_up_data_psd], group_labels_psd, 'Colors', 'rb', 'Symbol', '');
xtickangle(45);
title('PSD Correlation Comparison');
ylabel('Correlation Coefficient');
grid on;

% Coherence Plot
subplot(1,2,2);
hold on;
all_p_data_coh = [];
all_up_data_coh = [];
group_labels_coh = {};
for b = 1:length(metrics.coh)
    band = metrics.coh{b};
    p_data = arrayfun(@(x) x.coh.(band), paired_results);
    up_data = arrayfun(@(x) x.coh.(band), unpaired_results);

    all_p_data_coh = [all_p_data_coh, p_data];
    all_up_data_coh = [all_up_data_coh, up_data];
    group_labels_coh = [group_labels_coh, repmat({[band ' Paired']}, 1, length(p_data)), ...
                          repmat({[band ' Unpaired']}, 1, length(up_data))];

    [~, p] = ttest2(p_data, up_data);
    text(b*2-0.5, max([p_data, up_data])+0.02, sprintf('p=%.3f', p), ...
        'HorizontalAlignment', 'center');
end
boxplot([all_p_data_coh, all_up_data_coh], group_labels_coh, 'Colors', 'rb', 'Symbol', '');
xtickangle(45);
title('Coherence Comparison');
ylabel('Average Coherence');
grid on;


%% ============ SECTION 3.1: Signal-Level Comparisons ============
% Average PSD Spectra
if ~isempty(paired_results) && isfield(paired_results(1), 'psd_s_spectrum') && ~isempty(paired_results(1).psd_s_spectrum)
    figure('Name', 'Average PSD Spectrum Comparison', 'Position', [100, 100, 1000, 400]);
    subplot(1,2,1);
    hold on;

    paired_psd_s_cell = arrayfun(@(x) mean(x.psd_s_spectrum, 2), paired_results, 'UniformOutput', false);
    paired_psd_r_cell = arrayfun(@(x) mean(x.psd_r_spectrum, 2), paired_results, 'UniformOutput', false);

    unpaired_psd_s_cell = arrayfun(@(x) mean(x.psd_s_spectrum, 2), unpaired_results, 'UniformOutput', false);
    unpaired_psd_r_cell = arrayfun(@(x) mean(x.psd_r_spectrum, 2), unpaired_results, 'UniformOutput', false); 

    if ~isempty(paired_psd_s_cell)
        first_valid_idx_p = find(~cellfun('isempty', paired_psd_s_cell), 1);
        if isempty(first_valid_idx_p)
            warning('No valid PSD spectra found for Paired group to plot.');
        else
            num_freq_bins = size(paired_psd_s_cell{first_valid_idx_p}, 1);
            temp_paired_psd_s = zeros(num_freq_bins, length(paired_psd_s_cell));
            temp_paired_psd_r = zeros(num_freq_bins, length(paired_psd_r_cell));

            for i = 1:length(paired_psd_s_cell)
                if ~isempty(paired_psd_s_cell{i}) && size(paired_psd_s_cell{i}, 1) == num_freq_bins
                    temp_paired_psd_s(:, i) = paired_psd_s_cell{i};
                    temp_paired_psd_r(:, i) = paired_psd_r_cell{i};
                else
                    warning('Paired Results: PSD spectrum dimension mismatch or empty for pair %d. Skipped.', i);
                end
            end
            mean_paired_psd_s = mean(temp_paired_psd_s, 2);
            mean_paired_psd_r = mean(temp_paired_psd_r, 2);
            
            if ~isempty(paired_results(1).freq_psd) && (isequal(size(paired_results(1).freq_psd), size(mean_paired_psd_s)') || ...
               isequal(size(paired_results(1).freq_psd), size(mean_paired_psd_s)))
                plot(paired_results(1).freq_psd, mean_paired_psd_s, 'b', 'LineWidth', 1.5, 'DisplayName', 'Paired Sender PSD');
                plot(paired_results(1).freq_psd, mean_paired_psd_r, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Paired Receiver PSD');
            end
        end
    end

    if ~isempty(unpaired_psd_s_cell)
        first_valid_idx_up = find(~cellfun('isempty', unpaired_psd_s_cell), 1);
        if isempty(first_valid_idx_up)
            warning('No valid PSD spectra found for Unpaired group to plot.');
        else
            num_freq_bins = size(unpaired_psd_s_cell{first_valid_idx_up}, 1);
            temp_unpaired_psd_s = zeros(num_freq_bins, length(unpaired_psd_s_cell));
            temp_unpaired_psd_r = zeros(num_freq_bins, length(unpaired_psd_r_cell)); 

            for i = 1:length(unpaired_psd_s_cell)
                if ~isempty(unpaired_psd_s_cell{i}) && size(unpaired_psd_s_cell{i}, 1) == num_freq_bins
                    temp_unpaired_psd_s(:, i) = unpaired_psd_s_cell{i};
                    temp_unpaired_psd_r(:, i) = unpaired_psd_r_cell{i};
                else
                    warning('Unpaired Results: PSD spectrum dimension mismatch or empty for pair %d. Skipped.', i);
                end
            end
            mean_unpaired_psd_s = mean(temp_unpaired_psd_s, 2);
            mean_unpaired_psd_r = mean(temp_unpaired_psd_r, 2);

            if ~isempty(unpaired_results(1).freq_psd) && (isequal(size(unpaired_results(1).freq_psd), size(mean_unpaired_psd_s)') || ...
               isequal(size(unpaired_results(1).freq_psd), size(mean_unpaired_psd_s)))
                plot(unpaired_results(1).freq_psd, mean_unpaired_psd_s, 'r', 'LineWidth', 1.5, 'DisplayName', 'Unpaired Sender PSD');
                plot(unpaired_results(1).freq_psd, mean_unpaired_psd_r, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Unpaired Receiver PSD');
            end
        end
    end

    title('Average PSD Spectra (Channel-Averaged)');
    xlabel('Frequency (Hz)');
    ylabel('log_{10}(Power)');
    legend('Location', 'best');
    grid on;
    xlim([0.5 45]); % Relevant frequency range
end

% Average Coherence Spectra
if ~isempty(paired_results) && isfield(paired_results(1), 'coh_spectrum') && ~isempty(paired_results(1).coh_spectrum)
    subplot(1,2,2);
    hold on;
    
    paired_coh_specs_cell = arrayfun(@(x) x.coh_spectrum, paired_results, 'UniformOutput', false);
    unpaired_coh_specs_cell = arrayfun(@(x) x.coh_spectrum, unpaired_results, 'UniformOutput', false);

    if ~isempty(paired_coh_specs_cell)
        first_valid_idx_p_coh = find(~cellfun('isempty', paired_coh_specs_cell), 1);
        if isempty(first_valid_idx_p_coh)
             warning('No valid Coherence spectra found for Paired group to plot.');
        else
            num_coh_bins = size(paired_coh_specs_cell{first_valid_idx_p_coh}, 1);
            temp_paired_coh = zeros(num_coh_bins, length(paired_coh_specs_cell));
            for i = 1:length(paired_coh_specs_cell)
                if ~isempty(paired_coh_specs_cell{i}) && size(paired_coh_specs_cell{i}, 1) == num_coh_bins
                    temp_paired_coh(:, i) = paired_coh_specs_cell{i};
                else
                    warning('Paired Results: Coherence spectrum dimension mismatch or empty for pair %d. Skipped.', i);
                end
            end
            mean_paired_coh = mean(temp_paired_coh, 2);

            if ~isempty(paired_results(1).freq_coh) && (isequal(size(paired_results(1).freq_coh), size(mean_paired_coh)') || ...
               isequal(size(paired_results(1).freq_coh), size(mean_paired_coh)))
                plot(paired_results(1).freq_coh, mean_paired_coh, 'b', 'LineWidth', 1.5, 'DisplayName', 'Paired Coherence');
            end
        end
    end

    if ~isempty(unpaired_coh_specs_cell)
        first_valid_idx_up_coh = find(~cellfun('isempty', unpaired_coh_specs_cell), 1);
        if isempty(first_valid_idx_up_coh)
            warning('No valid Coherence spectra found for Unpaired group to plot.');
        else
            num_coh_bins = size(unpaired_coh_specs_cell{first_valid_idx_up_coh}, 1);
            temp_unpaired_coh = zeros(num_coh_bins, length(unpaired_coh_specs_cell));
            for i = 1:length(unpaired_coh_specs_cell)
                if ~isempty(unpaired_coh_specs_cell{i}) && size(unpaired_coh_specs_cell{i}, 1) == num_coh_bins
                    temp_unpaired_coh(:, i) = unpaired_coh_specs_cell{i};
                else
                    warning('Unpaired Results: Coherence spectrum dimension mismatch or empty for pair %d. Skipped.', i);
                end
            end
            mean_unpaired_coh = mean(temp_unpaired_coh, 2);

            if ~isempty(unpaired_results(1).freq_coh) && (isequal(size(unpaired_results(1).freq_coh), size(mean_unpaired_coh)') || ...
               isequal(size(unpaired_results(1).freq_coh), size(mean_unpaired_coh)))
                plot(unpaired_results(1).freq_coh, mean_unpaired_coh, 'r', 'LineWidth', 1.5, 'DisplayName', 'Unpaired Coherence');
            end
        end
    end

    title('Average Coherence Spectra');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    legend('Location', 'best');
    grid on;
    xlim([0.5 45]); % Relevant frequency range
end


%% ============ SECTION 3.2: Time Series Signal Comparison ============

% Average Epoch Signals (Channel-Averaged)
if ~isempty(paired_results) && isfield(paired_results(1), 'mean_s_epoch') && ~isempty(paired_results(1).mean_s_epoch)
    figure('Name', 'Average Epoch Signal Comparison', 'Position', [100, 100, 800, 400]);
    hold on;

    all_paired_s_epochs = cat(1, paired_results.mean_s_epoch);
    all_paired_r_epochs = cat(1, paired_results.mean_r_epoch);
    mean_paired_s_epoch_total = mean(all_paired_s_epochs, 1);
    mean_paired_r_epoch_total = mean(all_paired_r_epochs, 1);

    plot(paired_results(1).epoch_time, mean_paired_s_epoch_total, 'b', 'LineWidth', 1.5, 'DisplayName', 'Paired Sender Avg. Epoch');
    plot(paired_results(1).epoch_time, mean_paired_r_epoch_total, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Paired Receiver Avg. Epoch');

    all_unpaired_s_epochs = cat(1, unpaired_results.mean_s_epoch);
    all_unpaired_r_epochs = cat(1, unpaired_results.mean_r_epoch);
    mean_unpaired_s_epoch_total = mean(all_unpaired_s_epochs, 1);
    mean_unpaired_r_epoch_total = mean(all_unpaired_r_epochs, 1);

    plot(unpaired_results(1).epoch_time, mean_unpaired_s_epoch_total, 'r', 'LineWidth', 1.5, 'DisplayName', 'Unpaired Sender Avg. Epoch');
    plot(unpaired_results(1).epoch_time, mean_unpaired_r_epoch_total, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Unpaired Receiver Avg. Epoch');

    title('Average Epoch Signals (Channel-Averaged)');
    xlabel('Time (s)');
    ylabel('EEG Amplitude (Z-score)');
    legend('Location', 'best');
    grid on;
end


%% ============ SECTION 4: PERMUTATION TEST & EFFECT SIZE ============
fprintf('\nPermutation Test (Alpha Band Coherence Difference)...\n');
n_perm = 5000;
alpha_coh_paired = arrayfun(@(x) x.coh.Alpha, paired_results);
alpha_coh_unpaired = arrayfun(@(x) x.coh.Alpha, unpaired_results);
obs_diff = mean(alpha_coh_paired) - mean(alpha_coh_unpaired);
all_data = [alpha_coh_paired, alpha_coh_unpaired];
perm_diffs = zeros(n_perm, 1);
for k = 1:n_perm
    shuffled = all_data(randperm(length(all_data)));
    perm_diffs(k) = mean(shuffled(1:length(alpha_coh_paired))) - ...
                   mean(shuffled(length(alpha_coh_paired)+1:end));
end
p_perm = (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1);
fprintf('Permutation p-value: %.4f\n', p_perm);

if ~isempty(alpha_coh_paired) && ~isempty(alpha_coh_unpaired)
    cohens_d_alpha_coh = cohensd(alpha_coh_paired, alpha_coh_unpaired);
    fprintf('Cohen''s d (Alpha Coherence): %.3f\n', cohens_d_alpha_coh);
else
    fprintf('Not enough data to calculate Cohen''s d for Alpha Coherence.\n');
end


figure;
histogram(perm_diffs, 50, 'Normalization', 'probability', ...
         'FaceColor', [0.7 0.7 0.7]);
hold on;
xline(obs_diff, 'r-', 'LineWidth', 2);
xline(-obs_diff, 'r--', 'LineWidth', 1.5);
title('Permutation Distribution (Alpha Band Coherence Difference)');
xlabel('Mean Difference (Paired - Unpaired)');
ylabel('Probability Density');
legend('Null Distribution', sprintf('Observed: %.3f', obs_diff), ...
       sprintf('Symmetric: %.3f', -obs_diff));
grid on;


%% ============ SECTION 5: SEGMENT-BASED MEASUREMENTS SUMMARY TABLE ============
fprintf('\n=================================================================\n');
fprintf('                Segment-Based Measurements Summary               \n');
fprintf('=================================================================\n');

% Prepare data for table
metrics_table = cell(length(Config.bnames)*2, 5); % Rows for each band x 2 (PSD, Coherence), 5 columns
metrics_table_row_idx = 1;

fprintf('%-18s %-12s %-12s %-15s %-15s\n', 'Metric / Band', 'Paired (Mean)', 'Paired (Std)', 'Unpaired (Mean)', 'Unpaired (Std)');
fprintf('-----------------------------------------------------------------\n');

% PSD Correlation
for b_idx = 1:length(Config.bnames)
    band_name = Config.bnames{b_idx};
    
    p_data_psd = arrayfun(@(x) x.psd_corr.(band_name), paired_results);
    up_data_psd = arrayfun(@(x) x.psd_corr.(band_name), unpaired_results);
    
    mean_p_psd = mean(p_data_psd);
    std_p_psd = std(p_data_psd);
    mean_up_psd = mean(up_data_psd);
    std_up_psd = std(up_data_psd);
    
    fprintf('PSD_%-12s %-12.3f %-12.3f %-15.3f %-15.3f\n', band_name, mean_p_psd, std_p_psd, mean_up_psd, std_up_psd);
end

% Coherence
for b_idx = 1:length(Config.bnames)
    band_name = Config.bnames{b_idx};
    
    p_data_coh = arrayfun(@(x) x.coh.(band_name), paired_results);
    up_data_coh = arrayfun(@(x) x.coh.(band_name), unpaired_results);
    
    mean_p_coh = mean(p_data_coh);
    std_p_coh = std(p_data_coh);
    mean_up_coh = mean(up_data_coh);
    std_up_coh = std(up_data_coh);
    
    fprintf('COH_%-12s %-12.3f %-12.3f %-15.3f %-15.3f\n', band_name, mean_p_coh, std_p_coh, mean_up_coh, std_up_coh);
end
fprintf('=================================================================\n');



%% ================== HELPER FUNCTIONS ==================
% Function to segment EEG data
function [segmented_s_eeg, segmented_r_eeg, num_segments] = segment_eeg_data(eeg_s_data, eeg_r_data, segment_len_samples)
    data_len = size(eeg_s_data, 1);
    num_channels = size(eeg_s_data, 2);
    
    num_segments = floor(data_len / segment_len_samples);
    
    segmented_s_eeg = zeros(segment_len_samples, num_channels, num_segments);
    segmented_r_eeg = zeros(segment_len_samples, num_channels, num_segments);
    
    for i = 1:num_segments
        start_idx = (i - 1) * segment_len_samples + 1;
        end_idx = i * segment_len_samples;
        
        segmented_s_eeg(:, :, i) = eeg_s_data(start_idx:end_idx, :);
        segmented_r_eeg(:, :, i) = eeg_r_data(start_idx:end_idx, :);
    end
end

% preprocess_eeg, analyze_psd_correlation, analyze_coherence, cohensd functions remain unchanged
function eeg_clean = preprocess_eeg(eeg_raw, fs)
    % 1. Z-score normalization
    eeg_z = zscore(eeg_raw);
      
    % 2. 50 Hz notch filter
    [b,a] = iirnotch(50/(fs/2), 50/(fs/2)/35);
    eeg_notched = filtfilt(b, a, eeg_z);
      
    % 3. Bandpass filter (0.5-45 Hz)
    [b,a] = butter(4, [0.5 45]/(fs/2), 'bandpass');
    eeg_clean = filtfilt(b, a, eeg_notched);
end

% analyze_psd_correlation function updated: Returns spectra, frequency axis
function [res, log_psd_s_all, log_psd_r_all, f_psd] = analyze_psd_correlation(eeg_s, eeg_r, labels, cfg)
    % Note: 'labels' is passed but might be largely '1's if coming from segmented data
    % So spectral analysis is run on the entire segment
    
    [psd_s, f_psd] = pwelch(eeg_s, cfg.fft.win, cfg.fft.noverlap, cfg.fft.nfft, cfg.fs);
    [psd_r, ~] = pwelch(eeg_r, cfg.fft.win, cfg.fft.noverlap, cfg.fft.nfft, cfg.fs);
    
    log_psd_s_all = log10(psd_s);
    log_psd_r_all = log10(psd_r);
    
    res = struct();
    for i = 1:length(cfg.bnames)
        band = cfg.bands.(cfg.bnames{i});
        idx = f_psd >= band(1) & f_psd <= band(2);
        
        mean_s_band_per_channel = mean(log_psd_s_all(idx, :), 1);
        mean_r_band_per_channel = mean(log_psd_r_all(idx, :), 1);
        
        corr_mat = corrcoef(mean_s_band_per_channel, mean_r_band_per_channel);
        if ~isnan(corr_mat(1,2))
            res.(cfg.bnames{i}) = corr_mat(1, 2);
        else
            res.(cfg.bnames{i}) = 0; % Handle NaN if correlation cannot be computed
        end
    end
end

% analyze_coherence function updated: Returns spectra, frequency axis
function [res, Cxy_all, f_coh] = analyze_coherence(eeg_s, eeg_r, labels, cfg)
    % Note: 'labels' is passed but might be largely '1's if coming from segmented data
    % So spectral analysis is run on the entire segment
    n_ch = size(eeg_s, 2);
    
    band_coh = zeros(length(cfg.bnames), n_ch);
    Cxy_ch_all = zeros(cfg.fft.nfft/2 + 1, n_ch);
    f_coh = [];

    for ch = 1:n_ch
        % Ensure channels have data for coherence computation
        if isempty(eeg_s(:, ch)) || isempty(eeg_r(:, ch)) || all(isnan(eeg_s(:, ch))) || all(isnan(eeg_r(:, ch)))
            warning('Empty or NaN channel data for coherence computation in channel %d. Skipping.', ch);
            continue; 
        end

        [Cxy, f_coh_current] = mscohere(eeg_s(:, ch), eeg_r(:, ch), ...
                                cfg.fft.win, cfg.fft.noverlap, cfg.fft.nfft, cfg.fs);
        Cxy_ch_all(:, ch) = Cxy;

        if isempty(f_coh)
            f_coh = f_coh_current;
        end

        for b = 1:length(cfg.bnames)
            band = cfg.bands.(cfg.bnames{b});
            idx = f_coh_current >= band(1) & f_coh_current <= band(2);
            band_coh(b, ch) = mean(Cxy(idx));
        end
    end
    
    res = struct();
    for b = 1:length(cfg.bnames)
        res.(cfg.bnames{b}) = mean(band_coh(b, :));
    end
    Cxy_all = mean(Cxy_ch_all, 2);
    
    % Handle cases where Cxy_all might be all zeros or NaNs if no valid channels
    if all(isnan(Cxy_all)) || all(Cxy_all == 0)
        warning('Cxy_all is all NaNs or zeros after coherence computation.');
    end
end

function d = cohensd(x, y)
    n1 = length(x);
    n2 = length(y);
    s1 = std(x);
    s2 = std(y);
    
    if (n1 + n2 - 2) == 0 || (s1 == 0 && s2 == 0)
        d = 0; % Or NaN, depending on desired behavior for no variance
        warning('Cannot compute Cohen''s d: pooled standard deviation is zero or sample size too small.');
        return;
    end
    
    pooled_std = sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2));
    
    if pooled_std == 0
        d = 0; % Or NaN, if pooled standard deviation is still zero
        warning('Cannot compute Cohen''s d: pooled standard deviation is zero after calculation.');
        return;
    end
    
    d = (mean(x) - mean(y)) / pooled_std;
end