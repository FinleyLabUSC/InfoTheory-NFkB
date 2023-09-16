% Includes a proper intialization procedure
rng(12121995)
sigma = 1.25;
sim_size = 2e3;
timepoints_init = (0:1:9000)';
timepoints_main = (0:1:9000)';
parameters_0 = [ 3, 500, 69, 42, 67, 134, 67, 134, 134, 67, 0.002, 1000, ...
               2.5e-4, 0, 1.5e-3, 0.0, 1.67e-3, 1, 0, 1, 0, 0.01, 0, 1.7e-2, ...
               0.01, 9.4, 0.1*0.12, 0.2*0.015, 2.04e-2, 4.556, 0.00075, 0.67, 8e-5, 1125, ...
               2e-4, 1.5, 1.38e-2, 1.67e-2, 2.5e-5, 4.1e-3 ];

TRAF_maxima_agg = zeros(sim_size,1);           
RIP_maxima_agg  = zeros(sim_size,1);           
TAK_maxima_agg  = zeros(sim_size,1);     
IKK_maxima_agg  = zeros(sim_size,1);           
NFkB_maxima_agg  = zeros(sim_size,1);           
    
parameters_matrix = repmat(parameters_0, [sim_size,1]);

% Create a bimodal lognormal sample
positive_comp_weight = 0.5;
positive_comp = lognrnd( log(parameters_0(1)), 0.5, [ floor(sim_size*positive_comp_weight) ,1] );   % Sampling from a lognormal distribution 
negative_comp = lognrnd( -3, 0.5, [ sim_size-floor(sim_size*positive_comp_weight) ,1] );   % Sampling from a lognormal distribution 
Ant_concs = [positive_comp; negative_comp];
Ant_concs = Ant_concs(randperm(length(Ant_concs)));
% Finish creating a bimodal lognormal sample

parameters_matrix(:,2) = lognrnd( log(parameters_0(2)), sigma, [sim_size,1] );  % Randomizing only CAR init. conc. to measure information passed to TRAF2         
parfor idx = 1:sim_size
           
     param_vec = parameters_matrix(idx,:);
     param_vec(1) = 0;    
     [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
     species_init = species_out(end,:);
 
     param_vec(1) = Ant_concs(idx);
     species_init(1) = Ant_concs(idx);
     [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
 
     TRAF_maxima_agg(idx) = max(observables_out(:,1));

end
'DONE WITH TRAF'


parameters_matrix(:,3) = lognrnd( log(parameters_0(3)), sigma, [sim_size,1] );   % Additionally randomizing TARF2 to measure information passed to RIP1                  
parfor idx = 1:sim_size
           
     param_vec = parameters_matrix(idx,:);
     param_vec(1) = 0;    
     [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
     species_init = species_out(end,:);
 
     param_vec(1) = Ant_concs(idx);
     species_init(1) = Ant_concs(idx);
     [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
 
     RIP_maxima_agg(idx) = max(observables_out(:,2));

end
'DONE WITH RIP'


parameters_matrix(:,4) = lognrnd( log(parameters_0(4)), sigma, [sim_size,1] );  % Additionally randomizing RIP1 to measure information passed to TAB2/TAK1                  
parfor idx = 1:sim_size
           
     param_vec = parameters_matrix(idx,:);
     param_vec(1) = 0;    
     [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
     species_init = species_out(end,:);
 
     param_vec(1) = Ant_concs(idx);
     species_init(1) = Ant_concs(idx);
     [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
 
     TAK_maxima_agg(idx) = max(observables_out(:,5));

end
'DONE WITH TAB2/TAK1'


parameters_matrix(:,5) = lognrnd( log(parameters_0(5)), sigma, [sim_size,1] );
parameters_matrix(:,7) = lognrnd( log(parameters_0(7)), sigma, [sim_size,1] ); % Additionally randomizing TAB2 and TAK1 to measure information passed to NEMO/IKKbeta                   
parfor idx = 1:sim_size
           
     param_vec = parameters_matrix(idx,:);
     param_vec(1) = 0;    
     [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
     species_init = species_out(end,:);
 
     param_vec(1) = Ant_concs(idx);
     species_init(1) = Ant_concs(idx);
     [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
 
     IKK_maxima_agg(idx) = max(observables_out(:,8));

end
'DONE WITH NEMO/IKKbeta'


parameters_matrix(:,6) = lognrnd( log(parameters_0(6)), sigma, [sim_size,1] );
parameters_matrix(:,8) = lognrnd( log(parameters_0(8)), sigma, [sim_size,1] ); % Additionally randomizing NEMO and IKKbeta to measure information passed to NFkB                         
parfor idx = 1:sim_size
           
     param_vec = parameters_matrix(idx,:);
     param_vec(1) = 0;    
     [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
     species_init = species_out(end,:);
 
     param_vec(1) = Ant_concs(idx);
     species_init(1) = Ant_concs(idx);
     [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
 
     NFkB_maxima_agg(idx) = max(observables_out(:,end));

end
'DONE WITH NFKB'

parameters_matrix(:,9) = lognrnd( log(parameters_0(9)), sigma, [sim_size,1] );
parameters_matrix(:,10) = lognrnd( log(parameters_0(10)), sigma, [sim_size,1] ); % Additionally randomizing NFkB and IkBa to measure information carried by NFkB                         
parfor idx = 1:sim_size
           
     param_vec = parameters_matrix(idx,:);
     param_vec(1) = 0;    
     [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
     species_init = species_out(end,:);
     NFkB_maxima_agg_baseline(idx) = observables_out_init(end,end);
 
     param_vec(1) = Ant_concs(idx);
     species_init(1) = Ant_concs(idx);
     [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
 
     NFkB_maxima_agg_extend(idx) = max(observables_out(:,end));

end
'DONE WITH NFKB_extend'

save("mutinfo_125.mat", "Ant_concs", "TRAF_maxima_agg","RIP_maxima_agg",...
        "TAK_maxima_agg","IKK_maxima_agg","NFkB_maxima_agg","NFkB_maxima_agg_extend", ...
                "NFkB_maxima_agg_baseline")

