filename_list = ["mutinfo_025", "mutinfo_05", "mutinfo_075", "mutinfo_1", ...
                    "mutinfo_125", "mutinfo_15", "mutinfo_175", "mutinfo_2"];

for idx=1:8

    filename = filename_list(idx);
    load(filename)

    
    mutinfo(idx,1) = mi_cont_cont(Ant_concs, TRAF_maxima_agg, 5);
    mutinfo(idx,2) = mi_cont_cont(Ant_concs, RIP_maxima_agg,  5);
    mutinfo(idx,3) = mi_cont_cont(Ant_concs, TAK_maxima_agg,  5);
    mutinfo(idx,4) = mi_cont_cont(Ant_concs, IKK_maxima_agg,  5);
    mutinfo(idx,5) = mi_cont_cont(Ant_concs, NFkB_maxima_agg, 5);
    mutinfo(idx,6) = mi_cont_cont(Ant_concs, NFkB_maxima_agg_extend./NFkB_maxima_agg_baseline, 5);
    
end

figure
bar(mutinfo/log(2)) % Divide by log(2) to convert MI to a base 2 estimate
title("Mutual Information")
xticks([1 2 3 4 5 6 7 8])
xticklabels(["0.25", "0.50", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0"])
xlabel("Noise Level")
ylim([0 8.9])
ylabel("Mutual Information, bits")
set(gca,'FontName','Arial','fontsize',32, 'fontweight', 'bold')
