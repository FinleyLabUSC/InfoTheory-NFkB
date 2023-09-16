%% Load all the data from the simulations

filenames{1} = ["results_20_025", "results_20_05", "results_20_075",...
                 "results_20_1", "results_20_125", "results_20_15", ...
                 "results_20_175", "results_20_2", "results_20_5",];

filenames{2} = ["results_01_025", "results_01_05", "results_01_075",...
                 "results_01_1", "results_01_125", "results_01_15", ...
                 "results_01_175", "results_01_2", "results_01_5",];

filenames{3} = ["results_1_025", "results_1_05", "results_1_075",...
                   "results_1_1", "results_1_125", "results_1_15", ...
                   "results_1_175", "results_1_2", "results_1_5",];

filenames{4} = ["results_065_025", "results_065_05", "results_065_075",...
                   "results_065_1", "results_065_125", "results_065_15", ...
                   "results_065_175", "results_065_2", "results_065_5",];

filenames{5} = ["results_146_025", "results_146_05", "results_146_075",...
                   "results_146_1", "results_146_125", "results_146_15", ...
                   "results_146_175", "results_146_2", "results_146_5",];

filenames{6} = ["results_05_025", "results_05_05", "results_05_075",...
                   "results_05_1", "results_05_125", "results_05_15", ...
                   "results_05_175", "results_05_2", "results_05_5",];

filenames{7} = ["results_18_025", "results_18_05", "results_18_075",...
                   "results_18_1", "results_18_125", "results_18_15", ...
                   "results_18_175", "results_18_2", "results_18_5"];

filenames{8} = ["results_042_025", "results_042_05", "results_042_075",...
                   "results_042_1", "results_042_125", "results_042_15", ...
                   "results_042_175", "results_042_2", "results_042_5"];

filenames{9} = ["results_078_025", "results_078_05", "results_078_075",...
                   "results_078_1", "results_078_125", "results_078_15", ...
                   "results_078_175", "results_078_2", "results_078_5"];

filenames{10} = ["results_127_025", "results_127_05", "results_127_075",...
                   "results_127_1", "results_127_125", "results_127_15", ...
                   "results_127_175", "results_127_2", "results_127_5"];

filenames{11} = ["results_208_025", "results_208_05", "results_208_075",...
                   "results_208_1", "results_208_125", "results_208_15", ...
                   "results_208_175", "results_208_2", "results_208_5"];

filenames{12} = ["results_036_025", "results_036_05", "results_036_075",...
                   "results_036_1", "results_036_125", "results_036_15", ...
                   "results_036_175", "results_036_2", "results_036_5"];

filenames{13} = ["results_239_025", "results_239_05", "results_239_075",...
                   "results_239_1", "results_239_125", "results_239_15", ...
                   "results_239_175", "results_239_2", "results_239_5"];

filenames{14} = ["results_023_025", "results_023_05", "results_023_075",...
                   "results_023_1", "results_023_125", "results_023_15", ...
                   "results_023_175", "results_023_2", "results_023_5"];

filenames{15} = ["results_059_025", "results_059_05", "results_059_075",...
                   "results_059_1", "results_059_125", "results_059_15", ...
                   "results_059_175", "results_059_2", "results_059_5"];

filenames{16} = ["results_157_025", "results_157_05", "results_157_075",...
                   "results_157_1", "results_157_125", "results_157_15", ...
                   "results_157_175", "results_157_2", "results_157_5"];

filenames{17} = ["results_339_025", "results_339_05", "results_339_075",...
                   "results_339_1", "results_339_125", "results_339_15", ...
                   "results_339_175", "results_339_2", "results_339_5"];

filenames{18} = ["results_088_025", "results_088_05", "results_088_075",...
                   "results_088_1", "results_088_125", "results_088_15", ...
                   "results_088_175", "results_088_2", "results_088_5"];

filenames{19} = ["results_111_025", "results_111_05", "results_111_075",...
                   "results_111_1", "results_111_125", "results_111_15", ...
                   "results_111_175", "results_111_2", "results_111_5"];

filenames{20} = ["results_789_025", "results_789_05", "results_789_075",...
                   "results_789_1", "results_789_125", "results_789_15", ...
                   "results_789_175", "results_789_2", "results_789_5"];

filenames{21} = ["results_0_025", "results_0_05", "results_0_075",...
                   "results_0_1", "results_0_125", "results_0_15", ...
                   "results_0_175", "results_0_2", "results_0_5"];


mutinfo = zeros(6,8);
noise_levels = ["0.25", "0.50", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"];

NFkB_dump = [];
Ant_concs = [];

for jdx=1:8
    
   filename = filenames{21}(jdx);
   load(filename)
   NFkB_base{jdx} = makesUnique(NFkB_maxima_agg);

end


for idx=1:20
   for jdx=1:8
    
    filename = filenames{idx}(jdx);
    load(filename)
    NFkB_maxima{idx,jdx} = makesUnique(NFkB_maxima_agg)./NFkB_base{jdx};

   end
end


%% Construct the figure for 2 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}];    
    Ant_concs = [20*ones(size(NFkB_maxima{1,jdx})); 0.1*ones(size(NFkB_maxima{2,jdx}))];
    mutinfo(1,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

%     figure
%     histogram(NFkB_maxima{2,jdx},0:2:70)
%     hold on
%     histogram(NFkB_maxima{1,jdx},0:2:70)
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
    Ant_concs = [20*ones(size(NFkB_maxima{1,jdx})); 0.1*ones(size(NFkB_maxima{2,jdx})); 1.0*ones(size(NFkB_maxima{3,jdx}))];
    mutinfo(2,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 4 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{4,jdx}; NFkB_maxima{5,jdx}];    
    Ant_concs = [20*ones(size(NFkB_maxima{1,jdx})); 0.1*ones(size(NFkB_maxima{2,jdx})); 0.65*ones(size(NFkB_maxima{4,jdx})); 1.46*ones(size(NFkB_maxima{5,jdx}))];
    mutinfo(3,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 5 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{3,jdx}; NFkB_maxima{6,jdx}; NFkB_maxima{7,jdx}];    
    Ant_concs = [20*ones(size(NFkB_maxima{1,jdx})); 0.1*ones(size(NFkB_maxima{2,jdx})); 1.0*ones(size(NFkB_maxima{3,jdx})); 0.5*ones(size(NFkB_maxima{6,jdx})); 1.80*ones(size(NFkB_maxima{7,jdx}))];
    mutinfo(4,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 6 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{8,jdx}; NFkB_maxima{9,jdx}; NFkB_maxima{10,jdx}; NFkB_maxima{11,jdx}];    
    Ant_concs = [20*ones(size(NFkB_maxima{1,jdx})); 0.1*ones(size(NFkB_maxima{2,jdx})); 0.42*ones(size(NFkB_maxima{8,jdx})); 0.78*ones(size(NFkB_maxima{9,jdx})); 1.27*ones(size(NFkB_maxima{10,jdx})); 2.08*ones(size(NFkB_maxima{11,jdx}))];
    mutinfo(5,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate
        
end


%% Construct the figure for 7 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{1,jdx}; NFkB_maxima{2,jdx}; NFkB_maxima{4,jdx}; NFkB_maxima{5,jdx}; NFkB_maxima{12,jdx}; NFkB_maxima{13,jdx}; NFkB_maxima{3,jdx}];    
    Ant_concs = [20*ones(size(NFkB_maxima{1,jdx})); 0.1*ones(size(NFkB_maxima{2,jdx})); 0.65*ones(size(NFkB_maxima{4,jdx})); 1.46*ones(size(NFkB_maxima{5,jdx})); 0.36*ones(size(NFkB_maxima{12,jdx})); 2.39*ones(size(NFkB_maxima{13,jdx})); 1.0*ones(size(NFkB_maxima{3,jdx}))];
    mutinfo(6,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end



%% Construct the figure for 11 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{2,jdx}; NFkB_maxima{14,jdx}; NFkB_maxima{8,jdx}; NFkB_maxima{15,jdx}; NFkB_maxima{9,jdx}; NFkB_maxima{3,jdx}; NFkB_maxima{10,jdx}; NFkB_maxima{16,jdx}; NFkB_maxima{11,jdx}; NFkB_maxima{17,jdx}; NFkB_maxima{1,jdx}];    
    Ant_concs = [0.1*ones(size(NFkB_maxima{2,jdx})); 0.23*ones(size(NFkB_maxima{14,jdx})); 0.42*ones(size(NFkB_maxima{8,jdx})); 0.59*ones(size(NFkB_maxima{15,jdx})); 0.78*ones(size(NFkB_maxima{9,jdx})); 1.0*ones(size(NFkB_maxima{3,jdx})); 1.27*ones(size(NFkB_maxima{10,jdx})); 1.57*ones(size(NFkB_maxima{16,jdx})); 2.08*ones(size(NFkB_maxima{11,jdx})); 3.39*ones(size(NFkB_maxima{17,jdx})); 20*ones(size(NFkB_maxima{1,jdx}))];
    mutinfo(7,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end


%% Construct the figure for 20 inputs
for jdx=1:8
    
    NFkB_dump = [NFkB_maxima{2,jdx}; NFkB_maxima{14,jdx}; NFkB_maxima{8,jdx}; NFkB_maxima{15,jdx}; NFkB_maxima{9,jdx}; NFkB_maxima{3,jdx}; NFkB_maxima{10,jdx}; NFkB_maxima{16,jdx}; NFkB_maxima{11,jdx}; NFkB_maxima{17,jdx}; NFkB_maxima{1,jdx}; ...
                    NFkB_maxima{12,jdx}; NFkB_maxima{6,jdx}; NFkB_maxima{4,jdx}; NFkB_maxima{5,jdx}; NFkB_maxima{7,jdx}; NFkB_maxima{13,jdx}; NFkB_maxima{20,jdx}; NFkB_maxima{18,jdx}; NFkB_maxima{19,jdx}];    
    Ant_concs = [0.1*ones(size(NFkB_maxima{2,jdx})); 0.23*ones(size(NFkB_maxima{14,jdx})); 0.42*ones(size(NFkB_maxima{8,jdx})); 0.59*ones(size(NFkB_maxima{15,jdx})); 0.78*ones(size(NFkB_maxima{9,jdx})); 1.0*ones(size(NFkB_maxima{3,jdx})); ...
                    1.27*ones(size(NFkB_maxima{10,jdx})); 1.57*ones(size(NFkB_maxima{16,jdx})); 2.08*ones(size(NFkB_maxima{11,jdx})); 3.39*ones(size(NFkB_maxima{17,jdx})); 20*ones(size(NFkB_maxima{1,jdx})); ...
                       0.36*ones(size(NFkB_maxima{12,jdx})); 0.5*ones(size(NFkB_maxima{6,jdx})); 0.65*ones(size(NFkB_maxima{4,jdx})); 1.46*ones(size(NFkB_maxima{5,jdx})); 1.8*ones(size(NFkB_maxima{7,jdx})); 2.39*ones(size(NFkB_maxima{13,jdx})); 7.89*ones(size(NFkB_maxima{20,jdx}));...
                          0.88*ones(size(NFkB_maxima{18,jdx})); 1.11*ones(size(NFkB_maxima{19,jdx}));];
    mutinfo(8,jdx) = mi_discrete_cont(NFkB_dump, Ant_concs, 5)/log(2); % Divide by log(2) to convert MI to a base 2 estimate

end

figure
plot(mutinfo',"LineWidth",5)
title("Mutual Information")
xticks([1 2 3 4 5 6 7 8 9])
ylim([0 4])
xlim([1 8])
xticklabels(["0.25", "0.50", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"])
xlabel("Noise Level")
ylabel("Mutual Information, bits")
set(gca,'FontName','Menlo','fontsize',32,'fontweight','bold')
legend(["2 inputs", "3 inputs", "4 inputs", "5 inputs", "6 inputs", "7 inputs", ...
                            "11 inputs", "20 inputs"], 'NumColumns',2)

figure(2)
plot(max(mutinfo),'LineWidth', 5)
title("Channel Capacity")
xticks([1 2 3 4 5 6 7 8 9])
ylim([0 4])
xlim([1 8])
xticklabels(noise_levels)
xlabel("Noise Level")
ylabel("Channel Capacity, bits")
set(gca,'FontName','Menlo','fontsize',32,'fontweight','bold')
legend(["2 inputs", "3 inputs", "4 inputs", "5 inputs", "6 inputs", "7 inputs"])



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
