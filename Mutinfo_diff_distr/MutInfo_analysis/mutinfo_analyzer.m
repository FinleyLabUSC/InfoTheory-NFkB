filename_list = ["mutinfo_025", "mutinfo_05", "mutinfo_075",...
                    "mutinfo_1", "mutinfo_125" "mutinfo_15",...
                        "mutinfo_175", "mutinfo_2"];

folder_list = [ "../MutInfo_10_extend/", "../MutInfo_08_extend/",...
                "../MutInfo_05_extend/", "../MutInfo_02_extend/"];


for jdx=1:4
    for idx=1:8
    
        filename = folder_list(jdx)+filename_list(idx);
        load(filename) 

        mutinfo_abs(jdx,idx) = mi_cont_cont(Ant_concs, makesUnique(NFkB_maxima_agg), 5)/log(2);

    end
end

for jdx=1:4
    for idx=1:8
    
        filename = folder_list(jdx)+filename_list(idx);
        load(filename) 

        mutinfo_rel(jdx,idx) = mi_cont_cont(Ant_concs, makesUnique(NFkB_maxima_agg./NFkB_maxima_agg_baseline'), 5)/log(2);

    end
end

chan_cap_abs = [2.586695642718027,1.971919838604284,1.519361359447958,1.221176012096247,...
                    1.027862961905875,0.905167041429203,0.825359598102715,0.780237549624905];
chan_cap_rel = [2.586695642718027,1.971919838604284,1.519361359447958,1.220382529823758,...
                    1.02597142840782,0.890113251310335,0.758636111816182,0.68317389760535];

codeff_abs = mutinfo_abs./chan_cap_abs;
codeff_rel = mutinfo_rel./chan_cap_rel;

figure
title("MutInf, upstr., Abs.")
plot(mutinfo_abs', "LineWidth", 5)
xticks([1 2 3 4 5 6 7 8])
xticklabels(["0.25","0.50","0.75","1.0","1.25","1.5","1.75","2.0"])
xlabel("Noise Level")
yticks([0 0.5 1.0 1.5 2.0])
yticklabels(["0" "0.5" "1" "1.5" "2"])
ylim([0 2])
xlim([1 8])
ylabel("Mutual Information, bits")
set(gca,'FontName','Arial','fontsize',34,'fontweight','bold')

figure
title("MutInf, upstr., Rel.")
plot(mutinfo_rel', "LineWidth", 5)
xticks([1 2 3 4 5 6 7 8])
xticklabels(["0.25","0.50","0.75","1.0","1.25","1.5","1.75","2.0"])
xlabel("Noise Level")
yticks([0 0.5 1.0 1.5 2.0])
yticklabels(["0" "0.5" "1" "1.5" "2"])
ylim([0 2])
xlim([1 8])
ylabel("Mutual Information, bits")
set(gca,'FontName','Arial','fontsize',34,'fontweight','bold')

figure
title("CodeEff, upstr., Abs.")
plot(codeff_abs', "LineWidth", 5)
xticks([1 2 3 4 5 6 7 8])
xticklabels(["0.25","0.50","0.75","1.0","1.25","1.5","1.75","2.0"])
xlabel("Noise Level")
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2 ])
yticklabels(["0" "20%" "40%" "60%" "80%" "100%" "120%"])
ylim([0 1.3])
xlim([1 8])
ylabel("Coding Efficiency")
set(gca,'FontName','Arial','fontsize',34,'fontweight','bold')

figure
title("CodeEff, upstr., Rel.")
plot(codeff_rel', "LineWidth", 5)
xticks([1 2 3 4 5 6 7 8])
xticklabels(["0.25","0.50","0.75","1.0","1.25","1.5","1.75","2.0"])
xlabel("Noise Level")
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2 ])
yticklabels(["0" "20%" "40%" "60%" "80%" "100%" "120%"])
ylim([0 1.3])
xlim([1 8])
ylabel("Coding Efficiency")
set(gca,'FontName','Arial','fontsize',34,'fontweight','bold')


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

