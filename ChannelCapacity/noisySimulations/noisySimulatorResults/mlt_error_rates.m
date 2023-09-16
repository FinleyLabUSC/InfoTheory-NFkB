%% Load all the data for MLE estimation

filenames{1} = ["results_extend_70_025", "results_extend_70_05", "results_extend_70_075",...
                 "results_extend_70_1", "results_extend_70_125", "results_extend_70_15", ...
                 "results_extend_70_175", "results_extend_70_2", "results_extend_70_5",];

filenames{2} = ["results_extend_07_025", "results_extend_07_05", "results_extend_07_075",...
                 "results_extend_07_1", "results_extend_07_125", "results_extend_07_15", ...
                 "results_extend_07_175", "results_extend_07_2", "results_extend_07_5",];

filenames{3} = ["results_extend_0_025", "results_extend_0_05", "results_extend_0_075",...
                   "results_extend_0_1", "results_extend_0_125", "results_extend_0_15", ...
                   "results_extend_0_175", "results_extend_0_2", "results_extend_0_5"];



NFkB_dump = [];
Ant_concs = [];

for jdx=1:8
    
   filename = filenames{3}(jdx);
   load(filename)
   NFkB_base{jdx} = makesUnique(NFkB_maxima_agg_extend);

end

for idx=1:2
   for jdx=1:8
    
    filename = filenames{idx}(jdx);
    load(filename)
    NFkB_maxima{idx,jdx} = makesUnique(NFkB_maxima_agg_extend)./NFkB_base{jdx};

   end
end


%% Obtain MLT from the data for different noise levels

false_pos_agg = zeros(1,8);
false_neg_agg = zeros(1,8);
for idx=1:8
    
    figure(1)
    hist_low  = histogram(NFkB_maxima{2,idx}, 0:2:70);
    hist_low_probs = hist_low.Values;
    hold on 
    hist_high = histogram(NFkB_maxima{1,idx}, 0:2:70);
    hist_high_probs = hist_high.Values;

    thresh = find(hist_low_probs-hist_high_probs < 0,1)
    false_pos_agg(idx) = sum(hist_low_probs(thresh:end));
    false_neg_agg(idx) = sum(hist_high_probs(1:thresh-1));
    hold off

    ylim([0 1000])

end

figure(2)
plot(false_pos_agg/10,'LineWidth', 5, 'LineStyle',':')
hold on
plot(false_neg_agg/10, 'LineWidth', 5, 'LineStyle',':')
legend(["False Positive", "False Negative"])
title("False Estimates based on MLT")
xticks([1 2 3 4 5 6 7 8])
yticks([0 25 50 75 100])
ylim([0 100])
xlim([1 8])
xticklabels(["0.25", "0.50", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"])
xlabel("Noise Level")
ylabel("False Rate, %")
set(gca,'FontName','Arial','fontsize',32,'fontweight','bold')


%% Helper function to eliminate duplicates

function unique_vec = makesUnique(nonunique_vec)

    [unique_items, ~, ic] = unique(nonunique_vec);
    frequencies = accumarray(ic,1);
    non_unique = unique_items(frequencies>1);

    if length(non_unique)>0
        for jdx=1:length(non_unique)
            tracker = (nonunique_vec==non_unique(jdx)); 
            tracker_cumul = cumsum(tracker);
            tracker_cumul = tracker_cumul(randperm(length(tracker_cumul)));
            nonunique_vec(tracker) = nonunique_vec(tracker).*(1+tracker_cumul(tracker)*1e-15);
        end
    end

    unique_vec = nonunique_vec;

end
