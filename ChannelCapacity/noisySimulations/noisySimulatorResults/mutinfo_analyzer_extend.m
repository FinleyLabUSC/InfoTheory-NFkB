%% Load all the data from the simulations

filenames{1} = ["results_extend_70_025", "results_extend_70_05", "results_extend_70_075",...
                 "results_extend_70_1", "results_extend_70_125", "results_extend_70_15", ...
                 "results_extend_70_175", "results_extend_70_2", "results_extend_70_5",];

filenames{2} = ["results_extend_07_025", "results_extend_07_05", "results_extend_07_075",...
                 "results_extend_07_1", "results_extend_07_125", "results_extend_07_15", ...
                 "results_extend_07_175", "results_extend_07_2", "results_extend_07_5",];

filenames{3} = ["results_extend_11_025", "results_extend_11_05", "results_extend_11_075",...
                   "results_extend_11_1", "results_extend_11_125", "results_extend_11_15", ...
                   "results_extend_11_175", "results_extend_11_2", "results_extend_11_5",];

filenames{4} = ["results_extend_75_025", "results_extend_75_05", "results_extend_75_075",...
                   "results_extend_75_1", "results_extend_75_125", "results_extend_75_15", ...
                   "results_extend_75_175", "results_extend_75_2", "results_extend_75_5",];

filenames{5} = ["results_extend_145_025", "results_extend_145_05", "results_extend_145_075",...
                   "results_extend_145_1", "results_extend_145_125", "results_extend_145_15", ...
                   "results_extend_145_175", "results_extend_145_2", "results_extend_145_5",];

filenames{6} = ["results_extend_55_025", "results_extend_55_05", "results_extend_55_075",...
                   "results_extend_55_1", "results_extend_55_125", "results_extend_55_15", ...
                   "results_extend_55_175", "results_extend_55_2", "results_extend_55_5",];

filenames{7} = ["results_extend_165_025", "results_extend_165_05", "results_extend_165_075",...
                   "results_extend_165_1", "results_extend_165_125", "results_extend_165_15", ...
                   "results_extend_165_175", "results_extend_165_2", "results_extend_165_5"];

filenames{8} = ["results_extend_45_025", "results_extend_45_05", "results_extend_45_075",...
                   "results_extend_45_1", "results_extend_45_125", "results_extend_45_15", ...
                   "results_extend_45_175", "results_extend_45_2", "results_extend_45_5"];

filenames{9} = ["results_extend_9_025", "results_extend_9_05", "results_extend_9_075",...
                   "results_extend_9_1", "results_extend_9_125", "results_extend_9_15", ...
                   "results_extend_9_175", "results_extend_9_2", "results_extend_9_5"];

filenames{10} = ["results_extend_13_025", "results_extend_13_05", "results_extend_13_075",...
                   "results_extend_13_1", "results_extend_13_125", "results_extend_13_15", ...
                   "results_extend_13_175", "results_extend_13_2", "results_extend_13_5"];

filenames{11} = ["results_extend_18_025", "results_extend_18_05", "results_extend_18_075",...
                   "results_extend_18_1", "results_extend_18_125", "results_extend_18_15", ...
                   "results_extend_18_175", "results_extend_18_2", "results_extend_18_5"];

filenames{12} = ["results_extend_4_025", "results_extend_4_05", "results_extend_4_075",...
                   "results_extend_4_1", "results_extend_4_125", "results_extend_4_15", ...
                   "results_extend_4_175", "results_extend_4_2", "results_extend_4_5"];

filenames{13} = ["results_extend_19_025", "results_extend_19_05", "results_extend_19_075",...
                   "results_extend_19_1", "results_extend_19_125", "results_extend_19_15", ...
                   "results_extend_19_175", "results_extend_19_2", "results_extend_19_5"];

filenames{14} = ["results_extend_0_025", "results_extend_0_05", "results_extend_0_075",...
                   "results_extend_0_1", "results_extend_0_125", "results_extend_0_15", ...
                   "results_extend_0_175", "results_extend_0_2", "results_extend_0_5"];


mutinfo = zeros(6,8);
noise_levels = ["0.25", "0.50", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"];

NFkB_dump = [];
Ant_concs = [];

for jdx=1:8
    
   filename = filenames{14}(jdx);
   load(filename)
   NFkB_base{jdx} = makesUnique(NFkB_maxima_agg_extend);

end

for idx=1:13
   for jdx=1:8
    
    filename = filenames{idx}(jdx);
    load(filename)
    NFkB_maxima{idx,jdx} = makesUnique(NFkB_maxima_agg_extend)./NFkB_base{jdx};

   end
end


%% Construct the figure for 2 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}];    
    Ant_concs = [70*ones(size(NFkB_maxima{1,jdx})); 0.7*ones(size(NFkB_maxima{2,jdx}))];
    mutinfo(1,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

%     figure
%     histogram(NFkB_maxima{2,jdx},0:1:70)
%     hold on
%     histogram(NFkB_maxima{1,jdx},0:1:70)
%     legend(["Low Stim.", "High Stim."])
%     title("Dist. of NFkB_{max} for noise level: " + noise_levels(jdx))
%     xlabel("Concentration Values, molec./\mum^2")
%     ylabel("Number of cells")
%     xlim([0 70])
%     ylim([0 1700])
%     set(gca,'FontName','Menlo','fontsize',20)


    

end


%% Construct the figure for 3 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{3,jdx}];    
    Ant_concs = [70*ones(size(NFkB_maxima{1,jdx})); 0.7*ones(size(NFkB_maxima{2,jdx})); 11*ones(size(NFkB_maxima{3,jdx}))];
    mutinfo(2,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 4 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{4,jdx}; NFkB_maxima{5,jdx}];    
    Ant_concs = [70*ones(size(NFkB_maxima{1,jdx})); 0.7*ones(size(NFkB_maxima{2,jdx})); 7.5*ones(size(NFkB_maxima{4,jdx})); 14.5*ones(size(NFkB_maxima{5,jdx}))];
    mutinfo(3,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 5 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{3,jdx}; NFkB_maxima{6,jdx}; NFkB_maxima{7,jdx}];    
    Ant_concs = [70*ones(size(NFkB_maxima{1,jdx})); 0.7*ones(size(NFkB_maxima{2,jdx})); 11.0*ones(size(NFkB_maxima{3,jdx})); 5.5*ones(size(NFkB_maxima{6,jdx})); 16.5*ones(size(NFkB_maxima{7,jdx}))];
    mutinfo(4,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 6 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{8,jdx}; NFkB_maxima{9,jdx}; NFkB_maxima{10,jdx}; NFkB_maxima{11,jdx}];    
    Ant_concs = [70*ones(size(NFkB_maxima{1,jdx})); 0.7*ones(size(NFkB_maxima{2,jdx})); 4.5*ones(size(NFkB_maxima{8,jdx})); 9.0*ones(size(NFkB_maxima{9,jdx})); 13.0*ones(size(NFkB_maxima{10,jdx})); 18*ones(size(NFkB_maxima{11,jdx}))];
    mutinfo(5,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 7 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{4,jdx}; NFkB_maxima{5,jdx}; NFkB_maxima{12,jdx}; NFkB_maxima{13,jdx}; NFkB_maxima{3,jdx}];    
    Ant_concs = [70*ones(size(NFkB_maxima{1,jdx})); 0.7*ones(size(NFkB_maxima{2,jdx})); 7.5*ones(size(NFkB_maxima{4,jdx})); 14.5*ones(size(NFkB_maxima{5,jdx})); 4.0*ones(size(NFkB_maxima{12,jdx})); 19.0*ones(size(NFkB_maxima{13,jdx})); 11.0*ones(size(NFkB_maxima{3,jdx}))];
    mutinfo(6,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end

figure
plot(mutinfo',"LineWidth",5)
title("Mutual Information")
xticks([1 2 3 4 5 6 7 8])
ylim([0 4])
xlim([1 8])
xticklabels(["0.25", "0.50", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"])
xlabel("Noise Level")
ylabel("Mutual Information, bits")
set(gca,'FontName','Menlo','fontsize',32, 'fontweight', 'bold')
legend(["2 inputs", "3 inputs", "4 inputs", "5 inputs", "6 inputs", "7 inputs", "8 inputs", "20 inputs"])


figure
plot(max(mutinfo),'LineWidth', 5)
title("Channel Capacity")
xticks([1 2 3 4 5 6 7 8])
ylim([0 4])
xlim([1 8])
xticklabels(noise_levels)
xlabel("Noise Level")
ylabel("Channel Capacity, bits")
set(gca,'FontName','Menlo','fontsize',32, 'fontweight', 'bold')

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
