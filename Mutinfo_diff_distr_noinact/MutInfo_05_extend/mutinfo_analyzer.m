filename_list = ["mutinfo_025", "mutinfo_05", "mutinfo_075",...
                  "mutinfo_1", "mutinfo_125" "mutinfo_15",...
                   "mutinfo_175", "mutinfo_2"];
mutinfo = zeros(8,5);

for idx=1:8

    filename = filename_list(idx);
    load(filename)

    
    mutinfo(idx,1) = mi_cont_cont(Ant_concs, makesUnique(TRAF_maxima_agg), 5)/log(2);
    mutinfo(idx,2) = mi_cont_cont(Ant_concs, makesUnique(RIP_maxima_agg),  5)/log(2);
    mutinfo(idx,3) = mi_cont_cont(Ant_concs, makesUnique(TAK_maxima_agg),  5)/log(2);
    mutinfo(idx,4) = mi_cont_cont(Ant_concs, makesUnique(IKK_maxima_agg),  5)/log(2);
    mutinfo(idx,5) = mi_cont_cont(Ant_concs, makesUnique(NFkB_maxima_agg), 5)/log(2);
    mutinfo(idx,6) = mi_cont_cont(Ant_concs, makesUnique(NFkB_maxima_agg_extend'), 5)/log(2);
    
end


b = bar(mutinfo, 'FaceColor','flat') % Divide by log(2) to convert MI to a base 2 estimate
b(6).CData = [0.6350 0.0780 0.1840];
title("Mutual Information")
xticks([1 2 3 4 5 6 7 8])
xticklabels(["0.25","0.50","0.75","1.0","1.25","1.5","1.75","2.0"])
xlabel("Noise Level")
yticks(1:9)
ylim([0 8.9])
ylabel("Mutual Information, bits")
set(gca,'FontName','Menlo','fontsize',34,'fontweight','bold')

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

