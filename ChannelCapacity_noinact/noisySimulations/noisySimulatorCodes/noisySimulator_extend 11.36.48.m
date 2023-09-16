function noisySimulator_extend(ant_val, sigma)

    % Includes a proper intialization procedure
    rng(12121995)
    sim_size = 1e3;
    timepoints_init = (0:1:8000)';
    timepoints_main = (0:1:9000)';
    parameters_0 = [ ant_val, 500, 69, 42, 67, 134, 67, 134, 134, 67, 0.002, 1000, ...
                   2.5e-4, 0, 1.5e-3, 0.0, 1.67e-3, 1, 0, 1, 0, 0.01, 0, 1.7e-2, ...
                   0.01, 9.4, 0.1*0.12, 0.2*0.015, 2.04e-2, 4.556, 0.00075, 0.67, 8e-5, 1125, ...
                   2e-4, 1.5, 1.38e-2, 1.67e-2, 2.5e-5, 4.1e-3 ];  
    parameters_0(28) = 0;
               
    NFkB_maxima_agg_extend  = zeros(sim_size,1);           
        
    parameters_matrix = repmat(parameters_0, [sim_size,1]);
    
    parameters_matrix(:,2) = lognrnd( log(parameters_0(2)), sigma, [sim_size,1] );  % Randomizing only CAR init. conc. to measure information passed to TRAF2         
    parameters_matrix(:,3) = lognrnd( log(parameters_0(3)), sigma, [sim_size,1] );   % Additionally randomizing TARF2 to measure information passed to RIP1                     
    parameters_matrix(:,4) = lognrnd( log(parameters_0(4)), sigma, [sim_size,1] );  % Additionally randomizing RIP1 to measure information passed to TAB2/TAK1                  
    
    parameters_matrix(:,5) = lognrnd( log(parameters_0(5)), sigma, [sim_size,1] );
    parameters_matrix(:,7) = lognrnd( log(parameters_0(7)), sigma, [sim_size,1] ); % Additionally randomizing TAB2 and TAK1 to measure information passed to NEMO/IKKbeta                   
    
    parameters_matrix(:,6) = lognrnd( log(parameters_0(6)), sigma, [sim_size,1] );
    parameters_matrix(:,8) = lognrnd( log(parameters_0(8)), sigma, [sim_size,1] ); % Additionally randomizing NEMO and IKKbeta to measure information passed to NFkB             

    parameters_matrix(:,9) = lognrnd( log(parameters_0(9)), sigma, [sim_size,1] );
    parameters_matrix(:,10) = lognrnd( log(parameters_0(10)), sigma, [sim_size,1] ); % Additionally randomizing NEMO and IKKbeta to measure information in NFkB conc.                         
    parfor idx = 1:sim_size
               
         param_vec = parameters_matrix(idx,:);
         param_vec(1) = 0;    
         [~, ~, species_out, observables_out_init] = NFkB_model(timepoints_init, param_vec);
         species_init = species_out(end,:);
     
         param_vec(1) = ant_val;
         species_init(1) = ant_val;
         [~, tpts, ~, observables_out] = NFkB_model(timepoints_main, param_vec, species_init);
     
         NFkB_maxima_agg_extend(idx) = max(observables_out(:,end));
    
    end
    
    name_tag = strrep(num2str(ant_val), '.','') + "_" + strrep(num2str(sigma), '.','');
    save("results_extend_" + name_tag, "NFkB_maxima_agg_extend")

end

