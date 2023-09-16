function [err, timepoints, species_out, observables_out] = NFkB_model( timepoints, parameters, species_init, suppress_plot )
%NFKB_MODEL Integrate reaction network and plot observables.
%   Integrates the reaction network corresponding to the BioNetGen model
%   'NFkB_model' and then (optionally) plots the observable trajectories,
%   or species trajectories if no observables are defined. Trajectories are
%   generated using either default or user-defined parameters and initial
%   species values. Integration is performed by the MATLAB stiff solver
%   'ode15s'. NFKB_MODEL returns an error value, a vector of timepoints,
%   species trajectories, and observable trajectories.
%   
%   [err, timepoints, species_out, observables_out]
%        = NFkB_model( timepoints, species_init, parameters, suppress_plot )
%
%   INPUTS:
%   -------
%   species_init    : row vector of 83 initial species populations.
%   timepoints      : column vector of time points returned by integrator.
%   parameters      : row vector of 40 model parameters.
%   suppress_plot   : 0 if a plot is desired (default), 1 if plot is suppressed.
%
%   Note: to specify default value for an input argument, pass the empty array.
%
%   OUTPUTS:
%   --------
%   err             : 0 if the integrator exits without error, non-zero otherwise.
%   timepoints      : a row vector of timepoints returned by the integrator.
%   species_out     : array of species population trajectories
%                        (columns correspond to species, rows correspond to time).
%   observables_out : array of observable trajectories
%                        (columns correspond to observables, rows correspond to time).
%
%   QUESTIONS about the BNG Mfile generator?  Email justinshogg@gmail.com



%% Process input arguments

% define any missing arguments
if ( nargin < 1 )
    timepoints = [];
end

if ( nargin < 3 )
    species_init = [];
end

if ( nargin < 2 )
    parameters = [];
end

if ( nargin < 4 )
    suppress_plot = 1;
end


% initialize outputs (to avoid error msgs if script terminates early
err = 0;
species_out     = [];
observables_out = [];


% setup default parameters, if necessary
if ( isempty(parameters) )
   parameters = [ 1, 500, 6.9, 4.2, 6.7, 13.4, 6.7, 13.4, 13.4, 6.7, 1.7e-3, 1000, 2.5e-3, 0, 1.5e-2, 0, 1.67e-2, 1, 23.45, 1, 1072, 0.1, 4355, 2.5e-1, 0.1, 0.94, 0.25, 0.015, 2.04e-2, 0.4556, 0.0075, 0.067, 8e-5, 1125, 2e-4, 1.5, 1.38e-2, 1.67e-2, 2.5e-4, 4.1e-3 ];
end
% check that parameters has proper dimensions
if (  size(parameters,1) ~= 1  ||  size(parameters,2) ~= 40  )
    fprintf( 1, 'Error: size of parameter argument is invalid! Correct size = [1 40].\n' );
    err = 1;
    return;
end

% setup default initial values, if necessary
if ( isempty(species_init) )
   species_init = initialize_species( parameters );
end
% check that species_init has proper dimensions
if (  size(species_init,1) ~= 1  ||  size(species_init,2) ~= 83  )
    fprintf( 1, 'Error: size of species_init argument is invalid! Correct size = [1 83].\n' );
    err = 1;
    return;
end

% setup default timepoints, if necessary
if ( isempty(timepoints) )
   timepoints = linspace(0,10,20+1)';
end
% check that timepoints has proper dimensions
if (  size(timepoints,1) < 2  ||  size(timepoints,2) ~= 1  )
    fprintf( 1, 'Error: size of timepoints argument is invalid! Correct size = [t 1], t>1.\n' );
    err = 1;
    return;
end

% setup default suppress_plot, if necessary
if ( isempty(suppress_plot) )
   suppress_plot = 0;
end
% check that suppress_plot has proper dimensions
if ( size(suppress_plot,1) ~= 1  ||  size(suppress_plot,2) ~= 1 )
    fprintf( 1, 'Error: suppress_plots argument should be a scalar!\n' );
    err = 1;
    return;
end

% define parameter labels (this is for the user's reference!)
param_labels = { 'Ant_I', 'CAR_CD137_I', 'TRAF2_I', 'RIP1_I', 'TAB2_I', 'NEMO_I', 'TAK1_I', 'IKKb_I', 'IKBa_C_I', 'NFkB_C_I', 'k_ant_f', 'k_ant_d', 'k_car_cd137_traftri_f', 'k_car_cd137_traftri_d', 'k_traftri_rip1_f', 'k_traftri_rip1_d', 'k_rip1_ubq_f', 'k_rip1_tab2_f', 'k_rip1_tab2_d', 'k_rip1_nemo_f', 'k_rip1_nemo_d', 'k_tab2_tak1_f', 'k_tab2_tak1_d', 'k_tak1_act_f', 'k_nemo_ikkb_f', 'k_nemo_ikkb_d', 'k_ikkb_tak1', 'k_ikkb_inh', 'kcat_ikkb_ikba', 'km_ikkb', 'k_ikba_nfkb_f', 'k_ikba_nfkb_d', 'k_nfkb_exp', 'k_nfkb_eq', 'k_ikba_exp', 'k_ikba_eq', 'k_ikba_nfkb_exp', 'k_ikba_deg', 'k_transcript', 'k_translate' };



%% Integrate Network Model
 
% calculate expressions
[expressions] = calc_expressions( parameters );

% set ODE integrator options
opts = odeset( 'RelTol',   1e-8,   ...
               'AbsTol',   0.0001,   ...
               'Stats',    'off',  ...
               'BDF',      'off',    ...
               'MaxOrder', 5   );


% define derivative function
rhs_fcn = @(t,y)( calc_species_deriv( t, y, expressions ) );

% simulate model system (stiff integrator)
try 
    [~, species_out] = ode15s( rhs_fcn, timepoints, species_init', opts );
    if(length(timepoints) ~= size(species_out,1))
        exception = MException('ODE15sError:MissingOutput','Not all timepoints output\n');
        throw(exception);
    end
catch
    err = 1;
    fprintf( 1, 'Error: some problem encountered while integrating ODE network!\n' );
    return;
end

% calculate observables
observables_out = zeros( length(timepoints), 16 );
for t = 1 : length(timepoints)
    observables_out(t,:) = calc_observables( species_out(t,:), expressions );
end


%% Plot Output, if desired

if ( ~suppress_plot )
    
    % define plot labels
    observable_labels = { 'TRAF2_bound_C', 'RIP1_Ubq_C', 'TAB2_b_C', 'TAK1_I_C',...
                          'TAK1_A_C', 'NEMO_b_C', 'IKKb_N_C', 'IKKb_A_C', 'IKKb_I_C',...
                          'IKBa_C_C', 'IKBa_N_C', 'IKBaN_NFkB_C', 'IKBaP_NFkB_C',...
                          'IKBaN_NFkB_N_C', 'NFkB_cyto_C', 'NFkB_nucl_C' };

    % construct figure
    plot(timepoints,observables_out);
    title('NFkB_model observables','fontSize',14,'Interpreter','none');
    axis([0 timepoints(end) 0 inf]);
    legend(observable_labels,'fontSize',10,'Interpreter','none');
    xlabel('time','fontSize',12,'Interpreter','none');
    ylabel('number or concentration','fontSize',12,'Interpreter','none');

end


%~~~~~~~~~~~~~~~~~~~~~%
% END of main script! %
%~~~~~~~~~~~~~~~~~~~~~%

% Define if function to allow nested if statements in user-defined functions
function [val] = if__fun (cond, valT, valF)
% IF__FUN Select between two possible return values depending on the boolean
% variable COND.
    if (cond)
        val = valT;
    else
        val = valF;
    end
end

% initialize species function
function [species_init] = initialize_species( params )

    species_init = zeros(1,83);
    species_init(1) = params(1);
    species_init(2) = params(2);
    species_init(3) = params(3);
    species_init(4) = params(4);
    species_init(5) = params(5);
    species_init(6) = params(6);
    species_init(7) = params(7);
    species_init(8) = params(8);
    species_init(9) = params(10);
    species_init(10) = params(9);
    species_init(11) = 0;
    species_init(12) = 0;
    species_init(13) = 0;
    species_init(14) = 0;
    species_init(15) = 0;
    species_init(16) = 0;
    species_init(17) = 0;
    species_init(18) = 0;
    species_init(19) = 0;
    species_init(20) = 0;
    species_init(21) = 0;
    species_init(22) = 0;
    species_init(23) = 0;
    species_init(24) = 0;
    species_init(25) = 0;
    species_init(26) = 0;
    species_init(27) = 0;
    species_init(28) = 0;
    species_init(29) = 0;
    species_init(30) = 0;
    species_init(31) = 0;
    species_init(32) = 0;
    species_init(33) = 0;
    species_init(34) = 0;
    species_init(35) = 0;
    species_init(36) = 0;
    species_init(37) = 0;
    species_init(38) = 0;
    species_init(39) = 0;
    species_init(40) = 0;
    species_init(41) = 0;
    species_init(42) = 0;
    species_init(43) = 0;
    species_init(44) = 0;
    species_init(45) = 0;
    species_init(46) = 0;
    species_init(47) = 0;
    species_init(48) = 0;
    species_init(49) = 0;
    species_init(50) = 0;
    species_init(51) = 0;
    species_init(52) = 0;
    species_init(53) = 0;
    species_init(54) = 0;
    species_init(55) = 0;
    species_init(56) = 0;
    species_init(57) = 0;
    species_init(58) = 0;
    species_init(59) = 0;
    species_init(60) = 0;
    species_init(61) = 0;
    species_init(62) = 0;
    species_init(63) = 0;
    species_init(64) = 0;
    species_init(65) = 0;
    species_init(66) = 0;
    species_init(67) = 0;
    species_init(68) = 0;
    species_init(69) = 0;
    species_init(70) = 0;
    species_init(71) = 0;
    species_init(72) = 0;
    species_init(73) = 0;
    species_init(74) = 0;
    species_init(75) = 0;
    species_init(76) = 0;
    species_init(77) = 0;
    species_init(78) = 0;
    species_init(79) = 0;
    species_init(80) = 0;
    species_init(81) = 0;
    species_init(82) = 0;
    species_init(83) = 0;

end


% user-defined functions
% function rateLaw10
function [val] = rateLaw10(expressions, observables)
    val = (expressions(29)/(expressions(30)+observables(12)));
end




% Calculate expressions
function [ expressions ] = calc_expressions ( parameters )

    expressions = zeros(1,53);
    expressions(1) = parameters(1);
    expressions(2) = parameters(2);
    expressions(3) = parameters(3);
    expressions(4) = parameters(4);
    expressions(5) = parameters(5);
    expressions(6) = parameters(6);
    expressions(7) = parameters(7);
    expressions(8) = parameters(8);
    expressions(9) = parameters(9);
    expressions(10) = parameters(10);
    expressions(11) = parameters(11);
    expressions(12) = parameters(12);
    expressions(13) = parameters(13);
    expressions(14) = parameters(14);
    expressions(15) = parameters(15);
    expressions(16) = parameters(16);
    expressions(17) = parameters(17);
    expressions(18) = parameters(18);
    expressions(19) = parameters(19);
    expressions(20) = parameters(20);
    expressions(21) = parameters(21);
    expressions(22) = parameters(22);
    expressions(23) = parameters(23);
    expressions(24) = parameters(24);
    expressions(25) = parameters(25);
    expressions(26) = parameters(26);
    expressions(27) = parameters(27);
    expressions(28) = parameters(28);
    expressions(29) = parameters(29);
    expressions(30) = parameters(30);
    expressions(31) = parameters(31);
    expressions(32) = parameters(32);
    expressions(33) = parameters(33);
    expressions(34) = parameters(34);
    expressions(35) = parameters(35);
    expressions(36) = parameters(36);
    expressions(37) = parameters(37);
    expressions(38) = parameters(38);
    expressions(39) = parameters(39);
    expressions(40) = parameters(40);
    expressions(41) = (expressions(12)*expressions(11));
    expressions(42) = (expressions(14)*expressions(13));
    expressions(43) = (expressions(16)*expressions(15));
    expressions(44) = (expressions(19)*expressions(18));
    expressions(45) = (expressions(19)*expressions(18));
    expressions(46) = (expressions(21)*expressions(20));
    expressions(47) = (expressions(21)*expressions(20));
    expressions(48) = (expressions(23)*expressions(22));
    expressions(49) = (expressions(26)*expressions(25));
    expressions(50) = (expressions(32)*expressions(31));
    expressions(51) = (expressions(32)*expressions(31));
    expressions(52) = (expressions(36)*expressions(35));
    expressions(53) = (expressions(34)*expressions(33));
   
end



% Calculate observables
function [ observables ] = calc_observables ( species, expressions )

    observables = zeros(1,16);
    observables(1) = species(15) +species(18) +species(19) +species(20) +species(21) +species(22) +species(23) +species(24) +species(25) +species(26) +species(27) +species(28) +species(29) +species(30) +species(31) +species(32) +species(33) +species(34) +species(35) +species(36) +species(37) +species(38) +species(39) +species(40) +species(41) +species(42) +species(43) +species(44) +species(45) +species(46) +species(47) +species(48) +species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(55) +species(56) +species(57) +species(58) +species(59) +species(60) +species(61) +species(62) +species(63) +species(64) +species(65) +species(66) +species(67) +species(68) +species(69) +species(71) +species(72) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78) +species(79) +species(80) +species(81) +species(82) +species(83);
    observables(2) = species(19) +species(20) +species(21) +species(22) +species(23) +species(24) +species(25) +species(26) +species(27) +species(28) +species(29) +species(30) +species(31) +species(32) +species(33) +species(34) +species(35) +species(36) +species(37) +species(38) +species(39) +species(40) +species(41) +species(42) +species(43) +species(44) +species(45) +species(46) +species(47) +species(48) +species(49) +species(50) +species(51) +species(52) +species(53) +species(54) +species(55) +species(56) +species(57) +species(58) +species(59) +species(60) +species(61) +species(62) +species(63) +species(64) +species(65) +species(66) +species(67) +species(68) +species(69) +species(71) +species(72) +species(73) +species(74) +species(75) +species(76) +species(77) +species(78) +species(79) +species(80) +species(81) +species(82) +species(83);
    observables(3) = species(20) +species(21) +2*species(24) +species(25) +species(26) +species(28) +species(29) +2*species(32) +species(33) +2*species(34) +species(35) +species(36) +species(38) +species(40) +species(41) +2*species(42) +2*species(43) +species(44) +species(45) +2*species(46) +species(47) +species(48) +species(52) +species(53) +2*species(56) +2*species(57) +species(58) +species(59) +species(60) +species(61) +species(66) +species(67) +species(71) +species(72) +2*species(73) +species(74) +species(75) +species(79) +species(80);
    observables(4) = species(28) +species(29) +species(32) +species(34) +species(36) +species(38) +2*species(46) +species(47) +species(48) +species(56) +species(57) +species(58) +species(59) +species(71) +species(72);
    observables(5) = species(40) +species(41) +species(42) +species(43) +species(44) +species(45) +species(56) +species(57) +species(60) +species(61) +2*species(73) +species(74) +species(75) +species(79) +species(80);
    observables(6) = species(22) +species(23) +species(25) +species(26) +2*species(27) +species(30) +species(31) +species(33) +species(35) +species(36) +2*species(37) +species(38) +2*species(39) +species(44) +species(45) +species(47) +species(48) +2*species(49) +species(50) +species(51) +species(52) +species(53) +2*species(54) +2*species(55) +species(58) +species(59) +species(60) +species(61) +2*species(62) +2*species(63) +species(64) +species(65) +species(66) +species(67) +2*species(68) +2*species(69) +species(71) +species(72) +species(74) +species(75) +2*species(76) +2*species(77) +2*species(78) +species(79) +species(80) +2*species(81) +2*species(82) +2*species(83);
    observables(7) = species(30) +species(31) +species(33) +species(35) +species(37) +species(39) +species(47) +species(48) +2*species(49) +species(60) +species(61) +species(62) +species(63) +species(76) +species(77);
    observables(8) = species(50) +species(51) +species(52) +species(53) +species(54) +species(55) +species(58) +species(59) +species(62) +species(63) +species(74) +species(75) +2*species(78) +species(81) +species(82);
    observables(9) = species(64) +species(65) +species(66) +species(67) +species(68) +species(69) +species(71) +species(72) +species(76) +species(77) +species(79) +species(80) +species(81) +species(82) +2*species(83);
    observables(10) = species(10);
    observables(11) = species(13);
    observables(12) = species(12);
    observables(13) = species(70);
    observables(14) = species(16);
    observables(15) = species(9);
    observables(16) = species(14);

end


% Calculate ratelaws
function [ ratelaws ] = calcratelaws ( species, expressions, observables )

    ratelaws = zeros(1,16);
    ratelaws(1) = expressions(11)*species(1)*species(2);
    ratelaws(2) = expressions(31)*species(10)*species(9);
    ratelaws(3) = (expressions(36)*expressions(35))*species(10);
    ratelaws(4) = (expressions(34)*expressions(33))*species(9);
    ratelaws(5) = (expressions(12)*expressions(11))*species(11);
    ratelaws(6) = expressions(13)*species(11)*species(3);
    ratelaws(7) = (expressions(32)*expressions(31))*species(12);
    ratelaws(8) = expressions(31)*species(13)*species(14);
    ratelaws(9) = expressions(35)*species(13);
    ratelaws(10) = expressions(33)*species(14);
    ratelaws(11) = 0.5*expressions(39)*species(14)*species(14);
    ratelaws(12) = (expressions(14)*expressions(13))*species(15);
    ratelaws(13) = expressions(15)*species(15)*species(4);
    ratelaws(14) = (expressions(32)*expressions(31))*species(16);
    ratelaws(15) = expressions(37)*species(16);
    ratelaws(16) = expressions(40)*species(17);
    ratelaws(17) = (expressions(16)*expressions(15))*species(18);
    ratelaws(18) = expressions(17)*species(18);
    ratelaws(19) = expressions(18)*species(19)*species(5);
    ratelaws(20) = expressions(18)*species(19)*species(5);
    ratelaws(21) = expressions(20)*species(19)*species(6);
    ratelaws(22) = expressions(20)*species(19)*species(6);
    ratelaws(23) = expressions(18)*species(21)*species(5);
    ratelaws(24) = expressions(18)*species(23)*species(5);
    ratelaws(25) = (expressions(19)*expressions(18))*species(20);
    ratelaws(26) = expressions(18)*species(20)*species(5);
    ratelaws(27) = expressions(18)*species(22)*species(5);
    ratelaws(28) = (expressions(19)*expressions(18))*species(21);
    ratelaws(29) = expressions(20)*species(21)*species(6);
    ratelaws(30) = expressions(20)*species(23)*species(6);
    ratelaws(31) = (expressions(21)*expressions(20))*species(22);
    ratelaws(32) = expressions(20)*species(20)*species(6);
    ratelaws(33) = expressions(20)*species(22)*species(6);
    ratelaws(34) = (expressions(21)*expressions(20))*species(23);
    ratelaws(35) = expressions(22)*species(20)*species(7);
    ratelaws(36) = expressions(22)*species(21)*species(7);
    ratelaws(37) = expressions(25)*species(22)*species(8);
    ratelaws(38) = expressions(25)*species(23)*species(8);
    ratelaws(39) = expressions(18)*species(29)*species(5);
    ratelaws(40) = expressions(18)*species(31)*species(5);
    ratelaws(41) = (expressions(19)*expressions(18))*species(24);
    ratelaws(42) = (expressions(19)*expressions(18))*species(25);
    ratelaws(43) = expressions(18)*species(28)*species(5);
    ratelaws(44) = expressions(18)*species(30)*species(5);
    ratelaws(45) = (expressions(19)*expressions(18))*species(24);
    ratelaws(46) = (expressions(19)*expressions(18))*species(26);
    ratelaws(47) = expressions(20)*species(29)*species(6);
    ratelaws(48) = expressions(20)*species(31)*species(6);
    ratelaws(49) = (expressions(21)*expressions(20))*species(26);
    ratelaws(50) = (expressions(21)*expressions(20))*species(27);
    ratelaws(51) = expressions(20)*species(28)*species(6);
    ratelaws(52) = expressions(20)*species(30)*species(6);
    ratelaws(53) = (expressions(21)*expressions(20))*species(25);
    ratelaws(54) = (expressions(21)*expressions(20))*species(27);
    ratelaws(55) = expressions(22)*species(24)*species(7);
    ratelaws(56) = expressions(22)*species(24)*species(7);
    ratelaws(57) = expressions(22)*species(25)*species(7);
    ratelaws(58) = expressions(22)*species(26)*species(7);
    ratelaws(59) = (expressions(23)*expressions(22))*species(28);
    ratelaws(60) = (expressions(23)*expressions(22))*species(29);
    ratelaws(61) = expressions(24)*species(28);
    ratelaws(62) = expressions(24)*species(29);
    ratelaws(63) = expressions(25)*species(25)*species(8);
    ratelaws(64) = expressions(25)*species(26)*species(8);
    ratelaws(65) = expressions(25)*species(27)*species(8);
    ratelaws(66) = expressions(25)*species(27)*species(8);
    ratelaws(67) = (expressions(26)*expressions(25))*species(30);
    ratelaws(68) = (expressions(26)*expressions(25))*species(31);
    ratelaws(69) = expressions(18)*species(41)*species(5);
    ratelaws(70) = (expressions(19)*expressions(18))*species(32);
    ratelaws(71) = (expressions(19)*expressions(18))*species(33);
    ratelaws(72) = expressions(18)*species(40)*species(5);
    ratelaws(73) = (expressions(19)*expressions(18))*species(34);
    ratelaws(74) = (expressions(19)*expressions(18))*species(35);
    ratelaws(75) = expressions(20)*species(41)*species(6);
    ratelaws(76) = (expressions(21)*expressions(20))*species(36);
    ratelaws(77) = (expressions(21)*expressions(20))*species(37);
    ratelaws(78) = expressions(20)*species(40)*species(6);
    ratelaws(79) = (expressions(21)*expressions(20))*species(38);
    ratelaws(80) = (expressions(21)*expressions(20))*species(39);
    ratelaws(81) = expressions(22)*species(32)*species(7);
    ratelaws(82) = expressions(22)*species(33)*species(7);
    ratelaws(83) = expressions(22)*species(34)*species(7);
    ratelaws(84) = expressions(22)*species(35)*species(7);
    ratelaws(85) = (expressions(23)*expressions(22))*species(32);
    ratelaws(86) = (expressions(23)*expressions(22))*species(34);
    ratelaws(87) = (expressions(23)*expressions(22))*species(36);
    ratelaws(88) = (expressions(23)*expressions(22))*species(38);
    ratelaws(89) = expressions(24)*species(32);
    ratelaws(90) = expressions(24)*species(34);
    ratelaws(91) = expressions(24)*species(36);
    ratelaws(92) = expressions(24)*species(38);
    ratelaws(93) = expressions(25)*species(36)*species(8);
    ratelaws(94) = expressions(25)*species(37)*species(8);
    ratelaws(95) = expressions(25)*species(38)*species(8);
    ratelaws(96) = expressions(25)*species(39)*species(8);
    ratelaws(97) = (expressions(26)*expressions(25))*species(33);
    ratelaws(98) = (expressions(26)*expressions(25))*species(35);
    ratelaws(99) = (expressions(26)*expressions(25))*species(37);
    ratelaws(100) = (expressions(26)*expressions(25))*species(39);
    ratelaws(101) = expressions(27)*species(30)*species(40);
    ratelaws(102) = expressions(27)*species(30)*species(41);
    ratelaws(103) = expressions(27)*species(31)*species(40);
    ratelaws(104) = expressions(27)*species(31)*species(41);
    ratelaws(105) = expressions(27)*species(33)*species(40);
    ratelaws(106) = expressions(27)*species(33)*species(41);
    ratelaws(107) = expressions(27)*species(35)*species(40);
    ratelaws(108) = expressions(27)*species(35)*species(41);
    ratelaws(109) = expressions(27)*species(37)*species(40);
    ratelaws(110) = expressions(27)*species(37)*species(41);
    ratelaws(111) = expressions(27)*species(39)*species(40);
    ratelaws(112) = expressions(27)*species(39)*species(41);
    ratelaws(113) = expressions(18)*species(51)*species(5);
    ratelaws(114) = (expressions(19)*expressions(18))*species(42);
    ratelaws(115) = (expressions(19)*expressions(18))*species(52);
    ratelaws(116) = expressions(18)*species(50)*species(5);
    ratelaws(117) = (expressions(19)*expressions(18))*species(43);
    ratelaws(118) = (expressions(19)*expressions(18))*species(53);
    ratelaws(119) = expressions(20)*species(51)*species(6);
    ratelaws(120) = (expressions(21)*expressions(20))*species(44);
    ratelaws(121) = (expressions(21)*expressions(20))*species(54);
    ratelaws(122) = expressions(20)*species(50)*species(6);
    ratelaws(123) = (expressions(21)*expressions(20))*species(45);
    ratelaws(124) = (expressions(21)*expressions(20))*species(55);
    ratelaws(125) = expressions(22)*species(42)*species(7);
    ratelaws(126) = expressions(22)*species(43)*species(7);
    ratelaws(127) = expressions(22)*species(52)*species(7);
    ratelaws(128) = expressions(22)*species(53)*species(7);
    ratelaws(129) = (expressions(23)*expressions(22))*species(46);
    ratelaws(130) = (expressions(23)*expressions(22))*species(46);
    ratelaws(131) = (expressions(23)*expressions(22))*species(47);
    ratelaws(132) = (expressions(23)*expressions(22))*species(48);
    ratelaws(133) = expressions(24)*species(46);
    ratelaws(134) = expressions(24)*species(46);
    ratelaws(135) = expressions(24)*species(47);
    ratelaws(136) = expressions(24)*species(48);
    ratelaws(137) = expressions(25)*species(44)*species(8);
    ratelaws(138) = expressions(25)*species(45)*species(8);
    ratelaws(139) = expressions(25)*species(54)*species(8);
    ratelaws(140) = expressions(25)*species(55)*species(8);
    ratelaws(141) = (expressions(26)*expressions(25))*species(47);
    ratelaws(142) = (expressions(26)*expressions(25))*species(48);
    ratelaws(143) = (expressions(26)*expressions(25))*species(49);
    ratelaws(144) = (expressions(26)*expressions(25))*species(49);
    ratelaws(145) = expressions(27)*species(30)*species(42);
    ratelaws(146) = expressions(27)*species(30)*species(43);
    ratelaws(147) = expressions(27)*species(30)*species(44);
    ratelaws(148) = expressions(27)*species(30)*species(45);
    ratelaws(149) = expressions(27)*species(31)*species(42);
    ratelaws(150) = expressions(27)*species(31)*species(43);
    ratelaws(151) = expressions(27)*species(31)*species(44);
    ratelaws(152) = expressions(27)*species(31)*species(45);
    ratelaws(153) = expressions(27)*species(33)*species(42);
    ratelaws(154) = expressions(27)*species(33)*species(43);
    ratelaws(155) = expressions(27)*species(33)*species(44);
    ratelaws(156) = expressions(27)*species(33)*species(45);
    ratelaws(157) = expressions(27)*species(35)*species(42);
    ratelaws(158) = expressions(27)*species(35)*species(43);
    ratelaws(159) = expressions(27)*species(35)*species(44);
    ratelaws(160) = expressions(27)*species(35)*species(45);
    ratelaws(161) = expressions(27)*species(37)*species(42);
    ratelaws(162) = expressions(27)*species(37)*species(43);
    ratelaws(163) = expressions(27)*species(37)*species(44);
    ratelaws(164) = expressions(27)*species(37)*species(45);
    ratelaws(165) = expressions(27)*species(39)*species(42);
    ratelaws(166) = expressions(27)*species(39)*species(43);
    ratelaws(167) = expressions(27)*species(39)*species(44);
    ratelaws(168) = expressions(27)*species(39)*species(45);
    ratelaws(169) = expressions(27)*species(47)*species(40);
    ratelaws(170) = expressions(27)*species(47)*species(41);
    ratelaws(171) = expressions(27)*species(47)*species(42);
    ratelaws(172) = expressions(27)*species(47)*species(43);
    ratelaws(173) = expressions(27)*species(47)*species(44);
    ratelaws(174) = expressions(27)*species(47)*species(45);
    ratelaws(175) = expressions(27)*species(48)*species(40);
    ratelaws(176) = expressions(27)*species(48)*species(41);
    ratelaws(177) = expressions(27)*species(48)*species(42);
    ratelaws(178) = expressions(27)*species(48)*species(43);
    ratelaws(179) = expressions(27)*species(48)*species(44);
    ratelaws(180) = expressions(27)*species(48)*species(45);
    ratelaws(181) = expressions(27)*species(49)*species(40);
    ratelaws(182) = expressions(27)*species(49)*species(41);
    ratelaws(183) = expressions(27)*species(49)*species(42);
    ratelaws(184) = expressions(27)*species(49)*species(43);
    ratelaws(185) = expressions(27)*species(49)*species(44);
    ratelaws(186) = expressions(27)*species(49)*species(45);
    ratelaws(187) = expressions(27)*species(49)*species(40);
    ratelaws(188) = expressions(27)*species(49)*species(41);
    ratelaws(189) = expressions(27)*species(49)*species(42);
    ratelaws(190) = expressions(27)*species(49)*species(43);
    ratelaws(191) = expressions(27)*species(49)*species(44);
    ratelaws(192) = expressions(27)*species(49)*species(45);
    ratelaws(193) = expressions(28)*species(50);
    ratelaws(194) = expressions(28)*species(51);
    ratelaws(195) = expressions(28)*species(52);
    ratelaws(196) = expressions(28)*species(53);
    ratelaws(197) = expressions(28)*species(54);
    ratelaws(198) = expressions(28)*species(55);
    ratelaws(199) = rateLaw10(expressions,observables)*species(12)*species(50);
    ratelaws(200) = rateLaw10(expressions,observables)*species(12)*species(51);
    ratelaws(201) = rateLaw10(expressions,observables)*species(12)*species(52);
    ratelaws(202) = rateLaw10(expressions,observables)*species(12)*species(53);
    ratelaws(203) = rateLaw10(expressions,observables)*species(12)*species(54);
    ratelaws(204) = rateLaw10(expressions,observables)*species(12)*species(55);
    ratelaws(205) = expressions(18)*species(65)*species(5);
    ratelaws(206) = (expressions(19)*expressions(18))*species(66);
    ratelaws(207) = expressions(18)*species(64)*species(5);
    ratelaws(208) = (expressions(19)*expressions(18))*species(67);
    ratelaws(209) = expressions(20)*species(65)*species(6);
    ratelaws(210) = (expressions(21)*expressions(20))*species(68);
    ratelaws(211) = expressions(20)*species(64)*species(6);
    ratelaws(212) = (expressions(21)*expressions(20))*species(69);
    ratelaws(213) = expressions(22)*species(66)*species(7);
    ratelaws(214) = expressions(22)*species(67)*species(7);
    ratelaws(215) = (expressions(23)*expressions(22))*species(56);
    ratelaws(216) = (expressions(23)*expressions(22))*species(57);
    ratelaws(217) = (expressions(23)*expressions(22))*species(58);
    ratelaws(218) = (expressions(23)*expressions(22))*species(59);
    ratelaws(219) = expressions(24)*species(56);
    ratelaws(220) = expressions(24)*species(57);
    ratelaws(221) = expressions(24)*species(58);
    ratelaws(222) = expressions(24)*species(59);
    ratelaws(223) = expressions(25)*species(68)*species(8);
    ratelaws(224) = expressions(25)*species(69)*species(8);
    ratelaws(225) = (expressions(26)*expressions(25))*species(60);
    ratelaws(226) = (expressions(26)*expressions(25))*species(61);
    ratelaws(227) = (expressions(26)*expressions(25))*species(62);
    ratelaws(228) = (expressions(26)*expressions(25))*species(63);
    ratelaws(229) = expressions(27)*species(30)*species(56);
    ratelaws(230) = expressions(27)*species(30)*species(57);
    ratelaws(231) = expressions(27)*species(30)*species(60);
    ratelaws(232) = expressions(27)*species(30)*species(61);
    ratelaws(233) = expressions(27)*species(31)*species(56);
    ratelaws(234) = expressions(27)*species(31)*species(57);
    ratelaws(235) = expressions(27)*species(31)*species(60);
    ratelaws(236) = expressions(27)*species(31)*species(61);
    ratelaws(237) = expressions(27)*species(33)*species(56);
    ratelaws(238) = expressions(27)*species(33)*species(57);
    ratelaws(239) = expressions(27)*species(33)*species(60);
    ratelaws(240) = expressions(27)*species(33)*species(61);
    ratelaws(241) = expressions(27)*species(35)*species(56);
    ratelaws(242) = expressions(27)*species(35)*species(57);
    ratelaws(243) = expressions(27)*species(35)*species(60);
    ratelaws(244) = expressions(27)*species(35)*species(61);
    ratelaws(245) = expressions(27)*species(37)*species(56);
    ratelaws(246) = expressions(27)*species(37)*species(57);
    ratelaws(247) = expressions(27)*species(37)*species(60);
    ratelaws(248) = expressions(27)*species(37)*species(61);
    ratelaws(249) = expressions(27)*species(39)*species(56);
    ratelaws(250) = expressions(27)*species(39)*species(57);
    ratelaws(251) = expressions(27)*species(39)*species(60);
    ratelaws(252) = expressions(27)*species(39)*species(61);
    ratelaws(253) = expressions(27)*species(47)*species(56);
    ratelaws(254) = expressions(27)*species(47)*species(57);
    ratelaws(255) = expressions(27)*species(47)*species(60);
    ratelaws(256) = expressions(27)*species(47)*species(61);
    ratelaws(257) = expressions(27)*species(48)*species(56);
    ratelaws(258) = expressions(27)*species(48)*species(57);
    ratelaws(259) = expressions(27)*species(48)*species(60);
    ratelaws(260) = expressions(27)*species(48)*species(61);
    ratelaws(261) = expressions(27)*species(49)*species(56);
    ratelaws(262) = expressions(27)*species(49)*species(57);
    ratelaws(263) = expressions(27)*species(49)*species(60);
    ratelaws(264) = expressions(27)*species(49)*species(61);
    ratelaws(265) = expressions(27)*species(49)*species(56);
    ratelaws(266) = expressions(27)*species(49)*species(57);
    ratelaws(267) = expressions(27)*species(49)*species(60);
    ratelaws(268) = expressions(27)*species(49)*species(61);
    ratelaws(269) = expressions(27)*species(60)*species(40);
    ratelaws(270) = expressions(27)*species(60)*species(41);
    ratelaws(271) = expressions(27)*species(60)*species(42);
    ratelaws(272) = expressions(27)*species(60)*species(43);
    ratelaws(273) = expressions(27)*species(60)*species(44);
    ratelaws(274) = expressions(27)*species(60)*species(45);
    ratelaws(275) = expressions(27)*species(60)*species(56);
    ratelaws(276) = expressions(27)*species(60)*species(57);
    ratelaws(277) = expressions(27)*species(60)*species(60);
    ratelaws(278) = expressions(27)*species(60)*species(61);
    ratelaws(279) = expressions(27)*species(61)*species(40);
    ratelaws(280) = expressions(27)*species(61)*species(41);
    ratelaws(281) = expressions(27)*species(61)*species(42);
    ratelaws(282) = expressions(27)*species(61)*species(43);
    ratelaws(283) = expressions(27)*species(61)*species(44);
    ratelaws(284) = expressions(27)*species(61)*species(45);
    ratelaws(285) = expressions(27)*species(61)*species(56);
    ratelaws(286) = expressions(27)*species(61)*species(57);
    ratelaws(287) = expressions(27)*species(61)*species(60);
    ratelaws(288) = expressions(27)*species(61)*species(61);
    ratelaws(289) = expressions(27)*species(62)*species(40);
    ratelaws(290) = expressions(27)*species(62)*species(41);
    ratelaws(291) = expressions(27)*species(62)*species(42);
    ratelaws(292) = expressions(27)*species(62)*species(43);
    ratelaws(293) = expressions(27)*species(62)*species(44);
    ratelaws(294) = expressions(27)*species(62)*species(45);
    ratelaws(295) = expressions(27)*species(62)*species(56);
    ratelaws(296) = expressions(27)*species(62)*species(57);
    ratelaws(297) = expressions(27)*species(62)*species(60);
    ratelaws(298) = expressions(27)*species(62)*species(61);
    ratelaws(299) = expressions(27)*species(63)*species(40);
    ratelaws(300) = expressions(27)*species(63)*species(41);
    ratelaws(301) = expressions(27)*species(63)*species(42);
    ratelaws(302) = expressions(27)*species(63)*species(43);
    ratelaws(303) = expressions(27)*species(63)*species(44);
    ratelaws(304) = expressions(27)*species(63)*species(45);
    ratelaws(305) = expressions(27)*species(63)*species(56);
    ratelaws(306) = expressions(27)*species(63)*species(57);
    ratelaws(307) = expressions(27)*species(63)*species(60);
    ratelaws(308) = expressions(27)*species(63)*species(61);
    ratelaws(309) = expressions(28)*species(58);
    ratelaws(310) = expressions(28)*species(59);
    ratelaws(311) = expressions(28)*species(62);
    ratelaws(312) = expressions(28)*species(63);
    ratelaws(313) = rateLaw10(expressions,observables)*species(12)*species(58);
    ratelaws(314) = rateLaw10(expressions,observables)*species(12)*species(59);
    ratelaws(315) = rateLaw10(expressions,observables)*species(12)*species(62);
    ratelaws(316) = rateLaw10(expressions,observables)*species(12)*species(63);
    ratelaws(317) = expressions(38)*species(70);
    ratelaws(318) = (expressions(23)*expressions(22))*species(71);
    ratelaws(319) = (expressions(23)*expressions(22))*species(72);
    ratelaws(320) = expressions(24)*species(71);
    ratelaws(321) = expressions(24)*species(72);
    ratelaws(322) = (expressions(26)*expressions(25))*species(76);
    ratelaws(323) = (expressions(26)*expressions(25))*species(77);
    ratelaws(324) = expressions(27)*species(30)*species(73);
    ratelaws(325) = expressions(27)*species(30)*species(74);
    ratelaws(326) = expressions(27)*species(30)*species(75);
    ratelaws(327) = expressions(27)*species(31)*species(73);
    ratelaws(328) = expressions(27)*species(31)*species(74);
    ratelaws(329) = expressions(27)*species(31)*species(75);
    ratelaws(330) = expressions(27)*species(33)*species(73);
    ratelaws(331) = expressions(27)*species(33)*species(74);
    ratelaws(332) = expressions(27)*species(33)*species(75);
    ratelaws(333) = expressions(27)*species(35)*species(73);
    ratelaws(334) = expressions(27)*species(35)*species(74);
    ratelaws(335) = expressions(27)*species(35)*species(75);
    ratelaws(336) = expressions(27)*species(37)*species(73);
    ratelaws(337) = expressions(27)*species(37)*species(74);
    ratelaws(338) = expressions(27)*species(37)*species(75);
    ratelaws(339) = expressions(27)*species(39)*species(73);
    ratelaws(340) = expressions(27)*species(39)*species(74);
    ratelaws(341) = expressions(27)*species(39)*species(75);
    ratelaws(342) = expressions(27)*species(47)*species(73);
    ratelaws(343) = expressions(27)*species(47)*species(74);
    ratelaws(344) = expressions(27)*species(47)*species(75);
    ratelaws(345) = expressions(27)*species(48)*species(73);
    ratelaws(346) = expressions(27)*species(48)*species(74);
    ratelaws(347) = expressions(27)*species(48)*species(75);
    ratelaws(348) = expressions(27)*species(49)*species(73);
    ratelaws(349) = expressions(27)*species(49)*species(74);
    ratelaws(350) = expressions(27)*species(49)*species(75);
    ratelaws(351) = expressions(27)*species(49)*species(73);
    ratelaws(352) = expressions(27)*species(49)*species(74);
    ratelaws(353) = expressions(27)*species(49)*species(75);
    ratelaws(354) = expressions(27)*species(60)*species(73);
    ratelaws(355) = expressions(27)*species(60)*species(74);
    ratelaws(356) = expressions(27)*species(60)*species(75);
    ratelaws(357) = expressions(27)*species(61)*species(73);
    ratelaws(358) = expressions(27)*species(61)*species(74);
    ratelaws(359) = expressions(27)*species(61)*species(75);
    ratelaws(360) = expressions(27)*species(62)*species(73);
    ratelaws(361) = expressions(27)*species(62)*species(74);
    ratelaws(362) = expressions(27)*species(62)*species(75);
    ratelaws(363) = expressions(27)*species(63)*species(73);
    ratelaws(364) = expressions(27)*species(63)*species(74);
    ratelaws(365) = expressions(27)*species(63)*species(75);
    ratelaws(366) = expressions(27)*species(76)*species(40);
    ratelaws(367) = expressions(27)*species(76)*species(41);
    ratelaws(368) = expressions(27)*species(76)*species(42);
    ratelaws(369) = expressions(27)*species(76)*species(43);
    ratelaws(370) = expressions(27)*species(76)*species(44);
    ratelaws(371) = expressions(27)*species(76)*species(45);
    ratelaws(372) = expressions(27)*species(76)*species(56);
    ratelaws(373) = expressions(27)*species(76)*species(57);
    ratelaws(374) = expressions(27)*species(76)*species(60);
    ratelaws(375) = expressions(27)*species(76)*species(61);
    ratelaws(376) = expressions(27)*species(76)*species(73);
    ratelaws(377) = expressions(27)*species(76)*species(74);
    ratelaws(378) = expressions(27)*species(76)*species(75);
    ratelaws(379) = expressions(27)*species(77)*species(40);
    ratelaws(380) = expressions(27)*species(77)*species(41);
    ratelaws(381) = expressions(27)*species(77)*species(42);
    ratelaws(382) = expressions(27)*species(77)*species(43);
    ratelaws(383) = expressions(27)*species(77)*species(44);
    ratelaws(384) = expressions(27)*species(77)*species(45);
    ratelaws(385) = expressions(27)*species(77)*species(56);
    ratelaws(386) = expressions(27)*species(77)*species(57);
    ratelaws(387) = expressions(27)*species(77)*species(60);
    ratelaws(388) = expressions(27)*species(77)*species(61);
    ratelaws(389) = expressions(27)*species(77)*species(73);
    ratelaws(390) = expressions(27)*species(77)*species(74);
    ratelaws(391) = expressions(27)*species(77)*species(75);
    ratelaws(392) = expressions(28)*species(74);
    ratelaws(393) = expressions(28)*species(75);
    ratelaws(394) = expressions(28)*species(78);
    ratelaws(395) = expressions(28)*species(78);
    ratelaws(396) = rateLaw10(expressions,observables)*species(12)*species(74);
    ratelaws(397) = rateLaw10(expressions,observables)*species(12)*species(75);
    ratelaws(398) = rateLaw10(expressions,observables)*species(12)*species(78);
    ratelaws(399) = expressions(27)*species(30)*species(79);
    ratelaws(400) = expressions(27)*species(30)*species(80);
    ratelaws(401) = expressions(27)*species(31)*species(79);
    ratelaws(402) = expressions(27)*species(31)*species(80);
    ratelaws(403) = expressions(27)*species(33)*species(79);
    ratelaws(404) = expressions(27)*species(33)*species(80);
    ratelaws(405) = expressions(27)*species(35)*species(79);
    ratelaws(406) = expressions(27)*species(35)*species(80);
    ratelaws(407) = expressions(27)*species(37)*species(79);
    ratelaws(408) = expressions(27)*species(37)*species(80);
    ratelaws(409) = expressions(27)*species(39)*species(79);
    ratelaws(410) = expressions(27)*species(39)*species(80);
    ratelaws(411) = expressions(27)*species(47)*species(79);
    ratelaws(412) = expressions(27)*species(47)*species(80);
    ratelaws(413) = expressions(27)*species(48)*species(79);
    ratelaws(414) = expressions(27)*species(48)*species(80);
    ratelaws(415) = expressions(27)*species(49)*species(79);
    ratelaws(416) = expressions(27)*species(49)*species(80);
    ratelaws(417) = expressions(27)*species(49)*species(79);
    ratelaws(418) = expressions(27)*species(49)*species(80);
    ratelaws(419) = expressions(27)*species(60)*species(79);
    ratelaws(420) = expressions(27)*species(60)*species(80);
    ratelaws(421) = expressions(27)*species(61)*species(79);
    ratelaws(422) = expressions(27)*species(61)*species(80);
    ratelaws(423) = expressions(27)*species(62)*species(79);
    ratelaws(424) = expressions(27)*species(62)*species(80);
    ratelaws(425) = expressions(27)*species(63)*species(79);
    ratelaws(426) = expressions(27)*species(63)*species(80);
    ratelaws(427) = expressions(27)*species(76)*species(79);
    ratelaws(428) = expressions(27)*species(76)*species(80);
    ratelaws(429) = expressions(27)*species(77)*species(79);
    ratelaws(430) = expressions(27)*species(77)*species(80);
    ratelaws(431) = expressions(28)*species(81);
    ratelaws(432) = expressions(28)*species(82);
    ratelaws(433) = rateLaw10(expressions,observables)*species(12)*species(81);
    ratelaws(434) = rateLaw10(expressions,observables)*species(12)*species(82);

end

% Calculate species derivatives
function [ Dspecies ] = calc_species_deriv ( time, species, expressions )
    
    % initialize derivative vector
    Dspecies = zeros(83,1);
    
    % update observables
    [ observables ] = calc_observables( species, expressions );
    
    % update ratelaws
    [ ratelaws ] = calcratelaws( species, expressions, observables );
                        
    % calculate derivatives
    Dspecies(1) = -ratelaws(1) +ratelaws(5);
    Dspecies(2) = -ratelaws(1) +ratelaws(5);
    Dspecies(3) = -ratelaws(6) +ratelaws(12);
    Dspecies(4) = -ratelaws(13) +ratelaws(17);
    Dspecies(5) = -ratelaws(19) -ratelaws(20) -ratelaws(23) -ratelaws(24) +ratelaws(25) -ratelaws(26) -ratelaws(27) +ratelaws(28) -ratelaws(39) -ratelaws(40) +ratelaws(41) +ratelaws(42) -ratelaws(43) -ratelaws(44) +ratelaws(45) +ratelaws(46) -ratelaws(69) +ratelaws(70) +ratelaws(71) -ratelaws(72) +ratelaws(73) +ratelaws(74) -ratelaws(113) +ratelaws(114) +ratelaws(115) -ratelaws(116) +ratelaws(117) +ratelaws(118) -ratelaws(205) +ratelaws(206) -ratelaws(207) +ratelaws(208);
    Dspecies(6) = -ratelaws(21) -ratelaws(22) -ratelaws(29) -ratelaws(30) +ratelaws(31) -ratelaws(32) -ratelaws(33) +ratelaws(34) -ratelaws(47) -ratelaws(48) +ratelaws(49) +ratelaws(50) -ratelaws(51) -ratelaws(52) +ratelaws(53) +ratelaws(54) -ratelaws(75) +ratelaws(76) +ratelaws(77) -ratelaws(78) +ratelaws(79) +ratelaws(80) -ratelaws(119) +ratelaws(120) +ratelaws(121) -ratelaws(122) +ratelaws(123) +ratelaws(124) -ratelaws(209) +ratelaws(210) -ratelaws(211) +ratelaws(212);
    Dspecies(7) = -ratelaws(35) -ratelaws(36) -ratelaws(55) -ratelaws(56) -ratelaws(57) -ratelaws(58) +ratelaws(59) +ratelaws(60) -ratelaws(81) -ratelaws(82) -ratelaws(83) -ratelaws(84) +ratelaws(85) +ratelaws(86) +ratelaws(87) +ratelaws(88) -ratelaws(125) -ratelaws(126) -ratelaws(127) -ratelaws(128) +ratelaws(129) +ratelaws(130) +ratelaws(131) +ratelaws(132) -ratelaws(213) -ratelaws(214) +ratelaws(215) +ratelaws(216) +ratelaws(217) +ratelaws(218) +ratelaws(318) +ratelaws(319);
    Dspecies(8) = -ratelaws(37) -ratelaws(38) -ratelaws(63) -ratelaws(64) -ratelaws(65) -ratelaws(66) +ratelaws(67) +ratelaws(68) -ratelaws(93) -ratelaws(94) -ratelaws(95) -ratelaws(96) +ratelaws(97) +ratelaws(98) +ratelaws(99) +ratelaws(100) -ratelaws(137) -ratelaws(138) -ratelaws(139) -ratelaws(140) +ratelaws(141) +ratelaws(142) +ratelaws(143) +ratelaws(144) -ratelaws(223) -ratelaws(224) +ratelaws(225) +ratelaws(226) +ratelaws(227) +ratelaws(228) +ratelaws(322) +ratelaws(323);
    Dspecies(9) = -ratelaws(2) -ratelaws(4) +ratelaws(7) +ratelaws(10) +ratelaws(317);
    Dspecies(10) = -ratelaws(2) -ratelaws(3) +ratelaws(7) +ratelaws(9) +ratelaws(16);
    Dspecies(11) = ratelaws(1) -ratelaws(5) -ratelaws(6) +ratelaws(12);
    Dspecies(12) = ratelaws(2) -ratelaws(7) +ratelaws(15) -ratelaws(199) -ratelaws(200) -ratelaws(201) -ratelaws(202) -ratelaws(203) -ratelaws(204) -ratelaws(313) -ratelaws(314) -ratelaws(315) -ratelaws(316) -ratelaws(396) -ratelaws(397) -ratelaws(398) -ratelaws(433) -ratelaws(434);
    Dspecies(13) = ratelaws(3) -ratelaws(8) -ratelaws(9) +ratelaws(14);
    Dspecies(14) = ratelaws(4) -ratelaws(8) -ratelaws(10) +ratelaws(14);
    Dspecies(15) = ratelaws(6) -ratelaws(12) -ratelaws(13) +ratelaws(17);
    Dspecies(16) = ratelaws(8) -ratelaws(14) -ratelaws(15);
    Dspecies(17) = ratelaws(11) -ratelaws(16);
    Dspecies(18) = ratelaws(13) -ratelaws(17) -ratelaws(18);
    Dspecies(19) = ratelaws(18) -ratelaws(19) -ratelaws(20) -ratelaws(21) -ratelaws(22) +ratelaws(25) +ratelaws(28) +ratelaws(31) +ratelaws(34);
    Dspecies(20) = ratelaws(19) -ratelaws(25) -ratelaws(26) -ratelaws(32) -ratelaws(35) +ratelaws(45) +ratelaws(53) +ratelaws(59);
    Dspecies(21) = ratelaws(20) -ratelaws(23) -ratelaws(28) -ratelaws(29) -ratelaws(36) +ratelaws(41) +ratelaws(49) +ratelaws(60);
    Dspecies(22) = ratelaws(21) -ratelaws(27) -ratelaws(31) -ratelaws(33) -ratelaws(37) +ratelaws(46) +ratelaws(54) +ratelaws(67);
    Dspecies(23) = ratelaws(22) -ratelaws(24) -ratelaws(30) -ratelaws(34) -ratelaws(38) +ratelaws(42) +ratelaws(50) +ratelaws(68);
    Dspecies(24) = ratelaws(23) +ratelaws(26) -ratelaws(41) -ratelaws(45) -ratelaws(55) -ratelaws(56) +ratelaws(85) +ratelaws(86);
    Dspecies(25) = ratelaws(24) +ratelaws(32) -ratelaws(42) -ratelaws(53) -ratelaws(57) -ratelaws(63) +ratelaws(88) +ratelaws(97);
    Dspecies(26) = ratelaws(27) +ratelaws(29) -ratelaws(46) -ratelaws(49) -ratelaws(58) -ratelaws(64) +ratelaws(87) +ratelaws(98);
    Dspecies(27) = ratelaws(30) +ratelaws(33) -ratelaws(50) -ratelaws(54) -ratelaws(65) -ratelaws(66) +ratelaws(99) +ratelaws(100);
    Dspecies(28) = ratelaws(35) -ratelaws(43) -ratelaws(51) -ratelaws(59) -ratelaws(61) +ratelaws(73) +ratelaws(79);
    Dspecies(29) = ratelaws(36) -ratelaws(39) -ratelaws(47) -ratelaws(60) -ratelaws(62) +ratelaws(70) +ratelaws(76);
    Dspecies(30) = ratelaws(37) -ratelaws(44) -ratelaws(52) -ratelaws(67) +ratelaws(74) +ratelaws(80) -ratelaws(101) -ratelaws(102) -ratelaws(145) -ratelaws(146) -ratelaws(147) -ratelaws(148) -ratelaws(229) -ratelaws(230) -ratelaws(231) -ratelaws(232) -ratelaws(324) -ratelaws(325) -ratelaws(326) -ratelaws(399) -ratelaws(400);
    Dspecies(31) = ratelaws(38) -ratelaws(40) -ratelaws(48) -ratelaws(68) +ratelaws(71) +ratelaws(77) -ratelaws(103) -ratelaws(104) -ratelaws(149) -ratelaws(150) -ratelaws(151) -ratelaws(152) -ratelaws(233) -ratelaws(234) -ratelaws(235) -ratelaws(236) -ratelaws(327) -ratelaws(328) -ratelaws(329) -ratelaws(401) -ratelaws(402);
    Dspecies(32) = ratelaws(39) +ratelaws(55) -ratelaws(70) -ratelaws(81) -ratelaws(85) -ratelaws(89) +ratelaws(130);
    Dspecies(33) = ratelaws(40) +ratelaws(63) -ratelaws(71) -ratelaws(82) -ratelaws(97) -ratelaws(105) -ratelaws(106) +ratelaws(131) -ratelaws(153) -ratelaws(154) -ratelaws(155) -ratelaws(156) -ratelaws(237) -ratelaws(238) -ratelaws(239) -ratelaws(240) -ratelaws(330) -ratelaws(331) -ratelaws(332) -ratelaws(403) -ratelaws(404);
    Dspecies(34) = ratelaws(43) +ratelaws(56) -ratelaws(73) -ratelaws(83) -ratelaws(86) -ratelaws(90) +ratelaws(129);
    Dspecies(35) = ratelaws(44) +ratelaws(64) -ratelaws(74) -ratelaws(84) -ratelaws(98) -ratelaws(107) -ratelaws(108) +ratelaws(132) -ratelaws(157) -ratelaws(158) -ratelaws(159) -ratelaws(160) -ratelaws(241) -ratelaws(242) -ratelaws(243) -ratelaws(244) -ratelaws(333) -ratelaws(334) -ratelaws(335) -ratelaws(405) -ratelaws(406);
    Dspecies(36) = ratelaws(47) +ratelaws(58) -ratelaws(76) -ratelaws(87) -ratelaws(91) -ratelaws(93) +ratelaws(142);
    Dspecies(37) = ratelaws(48) +ratelaws(65) -ratelaws(77) -ratelaws(94) -ratelaws(99) -ratelaws(109) -ratelaws(110) +ratelaws(144) -ratelaws(161) -ratelaws(162) -ratelaws(163) -ratelaws(164) -ratelaws(245) -ratelaws(246) -ratelaws(247) -ratelaws(248) -ratelaws(336) -ratelaws(337) -ratelaws(338) -ratelaws(407) -ratelaws(408);
    Dspecies(38) = ratelaws(51) +ratelaws(57) -ratelaws(79) -ratelaws(88) -ratelaws(92) -ratelaws(95) +ratelaws(141);
    Dspecies(39) = ratelaws(52) +ratelaws(66) -ratelaws(80) -ratelaws(96) -ratelaws(100) -ratelaws(111) -ratelaws(112) +ratelaws(143) -ratelaws(165) -ratelaws(166) -ratelaws(167) -ratelaws(168) -ratelaws(249) -ratelaws(250) -ratelaws(251) -ratelaws(252) -ratelaws(339) -ratelaws(340) -ratelaws(341) -ratelaws(409) -ratelaws(410);
    Dspecies(40) = ratelaws(61) -ratelaws(72) -ratelaws(78) +ratelaws(117) +ratelaws(123);
    Dspecies(41) = ratelaws(62) -ratelaws(69) -ratelaws(75) +ratelaws(114) +ratelaws(120);
    Dspecies(42) = ratelaws(69) +ratelaws(89) -ratelaws(114) -ratelaws(125) +ratelaws(215);
    Dspecies(43) = ratelaws(72) +ratelaws(90) -ratelaws(117) -ratelaws(126) +ratelaws(216);
    Dspecies(44) = ratelaws(75) +ratelaws(91) -ratelaws(120) -ratelaws(137) +ratelaws(226);
    Dspecies(45) = ratelaws(78) +ratelaws(92) -ratelaws(123) -ratelaws(138) +ratelaws(225);
    Dspecies(46) = ratelaws(81) +ratelaws(83) -ratelaws(129) -ratelaws(130) -ratelaws(133) -ratelaws(134);
    Dspecies(47) = ratelaws(82) +ratelaws(95) -ratelaws(131) -ratelaws(135) -ratelaws(141) -ratelaws(169) -ratelaws(170) -ratelaws(171) -ratelaws(172) -ratelaws(173) -ratelaws(174) -ratelaws(253) -ratelaws(254) -ratelaws(255) -ratelaws(256) -ratelaws(342) -ratelaws(343) -ratelaws(344) -ratelaws(411) -ratelaws(412);
    Dspecies(48) = ratelaws(84) +ratelaws(93) -ratelaws(132) -ratelaws(136) -ratelaws(142) -ratelaws(175) -ratelaws(176) -ratelaws(177) -ratelaws(178) -ratelaws(179) -ratelaws(180) -ratelaws(257) -ratelaws(258) -ratelaws(259) -ratelaws(260) -ratelaws(345) -ratelaws(346) -ratelaws(347) -ratelaws(413) -ratelaws(414);
    Dspecies(49) = ratelaws(94) +ratelaws(96) -ratelaws(143) -ratelaws(144) -ratelaws(181) -ratelaws(182) -ratelaws(183) -ratelaws(184) -ratelaws(185) -ratelaws(186) -ratelaws(187) -ratelaws(188) -ratelaws(189) -ratelaws(190) -ratelaws(191) -ratelaws(192) -ratelaws(261) -ratelaws(262) -ratelaws(263) -ratelaws(264) -ratelaws(265) -ratelaws(266) -ratelaws(267) -ratelaws(268) -ratelaws(348) -ratelaws(349) -ratelaws(350) -ratelaws(351) -ratelaws(352) -ratelaws(353) -ratelaws(415) -ratelaws(416) -ratelaws(417) -ratelaws(418);
    Dspecies(50) = ratelaws(101) +ratelaws(102) -ratelaws(116) +ratelaws(118) -ratelaws(122) +ratelaws(124) +ratelaws(145) +ratelaws(146) +ratelaws(147) +ratelaws(148) -ratelaws(193) +ratelaws(229) +ratelaws(230) +ratelaws(231) +ratelaws(232) +ratelaws(324) +ratelaws(325) +ratelaws(326) +ratelaws(399) +ratelaws(400);
    Dspecies(51) = ratelaws(103) +ratelaws(104) -ratelaws(113) +ratelaws(115) -ratelaws(119) +ratelaws(121) +ratelaws(149) +ratelaws(150) +ratelaws(151) +ratelaws(152) -ratelaws(194) +ratelaws(233) +ratelaws(234) +ratelaws(235) +ratelaws(236) +ratelaws(327) +ratelaws(328) +ratelaws(329) +ratelaws(401) +ratelaws(402);
    Dspecies(52) = ratelaws(105) +ratelaws(106) +ratelaws(113) -ratelaws(115) -ratelaws(127) +ratelaws(153) +ratelaws(154) +ratelaws(155) +ratelaws(156) -ratelaws(195) +ratelaws(217) +ratelaws(237) +ratelaws(238) +ratelaws(239) +ratelaws(240) +ratelaws(330) +ratelaws(331) +ratelaws(332) +ratelaws(403) +ratelaws(404);
    Dspecies(53) = ratelaws(107) +ratelaws(108) +ratelaws(116) -ratelaws(118) -ratelaws(128) +ratelaws(157) +ratelaws(158) +ratelaws(159) +ratelaws(160) -ratelaws(196) +ratelaws(218) +ratelaws(241) +ratelaws(242) +ratelaws(243) +ratelaws(244) +ratelaws(333) +ratelaws(334) +ratelaws(335) +ratelaws(405) +ratelaws(406);
    Dspecies(54) = ratelaws(109) +ratelaws(110) +ratelaws(119) -ratelaws(121) -ratelaws(139) +ratelaws(161) +ratelaws(162) +ratelaws(163) +ratelaws(164) -ratelaws(197) +ratelaws(227) +ratelaws(245) +ratelaws(246) +ratelaws(247) +ratelaws(248) +ratelaws(336) +ratelaws(337) +ratelaws(338) +ratelaws(407) +ratelaws(408);
    Dspecies(55) = ratelaws(111) +ratelaws(112) +ratelaws(122) -ratelaws(124) -ratelaws(140) +ratelaws(165) +ratelaws(166) +ratelaws(167) +ratelaws(168) -ratelaws(198) +ratelaws(228) +ratelaws(249) +ratelaws(250) +ratelaws(251) +ratelaws(252) +ratelaws(339) +ratelaws(340) +ratelaws(341) +ratelaws(409) +ratelaws(410);
    Dspecies(56) = ratelaws(125) +ratelaws(133) -ratelaws(215) -ratelaws(219);
    Dspecies(57) = ratelaws(126) +ratelaws(134) -ratelaws(216) -ratelaws(220);
    Dspecies(58) = ratelaws(127) +ratelaws(169) +ratelaws(170) +ratelaws(171) +ratelaws(172) +ratelaws(173) +ratelaws(174) -ratelaws(217) -ratelaws(221) +ratelaws(253) +ratelaws(254) +ratelaws(255) +ratelaws(256) -ratelaws(309) +ratelaws(342) +ratelaws(343) +ratelaws(344) +ratelaws(411) +ratelaws(412);
    Dspecies(59) = ratelaws(128) +ratelaws(175) +ratelaws(176) +ratelaws(177) +ratelaws(178) +ratelaws(179) +ratelaws(180) -ratelaws(218) -ratelaws(222) +ratelaws(257) +ratelaws(258) +ratelaws(259) +ratelaws(260) -ratelaws(310) +ratelaws(345) +ratelaws(346) +ratelaws(347) +ratelaws(413) +ratelaws(414);
    Dspecies(60) = ratelaws(135) +ratelaws(138) -ratelaws(225) -ratelaws(269) -ratelaws(270) -ratelaws(271) -ratelaws(272) -ratelaws(273) -ratelaws(274) -ratelaws(275) -ratelaws(276) -ratelaws(277) -ratelaws(278) -ratelaws(354) -ratelaws(355) -ratelaws(356) -ratelaws(419) -ratelaws(420);
    Dspecies(61) = ratelaws(136) +ratelaws(137) -ratelaws(226) -ratelaws(279) -ratelaws(280) -ratelaws(281) -ratelaws(282) -ratelaws(283) -ratelaws(284) -ratelaws(285) -ratelaws(286) -ratelaws(287) -ratelaws(288) -ratelaws(357) -ratelaws(358) -ratelaws(359) -ratelaws(421) -ratelaws(422);
    Dspecies(62) = ratelaws(139) +ratelaws(181) +ratelaws(182) +ratelaws(183) +ratelaws(184) +ratelaws(185) +ratelaws(186) -ratelaws(227) +ratelaws(261) +ratelaws(262) +ratelaws(263) +ratelaws(264) -ratelaws(289) -ratelaws(290) -ratelaws(291) -ratelaws(292) -ratelaws(293) -ratelaws(294) -ratelaws(295) -ratelaws(296) -ratelaws(297) -ratelaws(298) -ratelaws(311) +ratelaws(348) +ratelaws(349) +ratelaws(350) -ratelaws(360) -ratelaws(361) -ratelaws(362) +ratelaws(415) +ratelaws(416) -ratelaws(423) -ratelaws(424);
    Dspecies(63) = ratelaws(140) +ratelaws(187) +ratelaws(188) +ratelaws(189) +ratelaws(190) +ratelaws(191) +ratelaws(192) -ratelaws(228) +ratelaws(265) +ratelaws(266) +ratelaws(267) +ratelaws(268) -ratelaws(299) -ratelaws(300) -ratelaws(301) -ratelaws(302) -ratelaws(303) -ratelaws(304) -ratelaws(305) -ratelaws(306) -ratelaws(307) -ratelaws(308) -ratelaws(312) +ratelaws(351) +ratelaws(352) +ratelaws(353) -ratelaws(363) -ratelaws(364) -ratelaws(365) +ratelaws(417) +ratelaws(418) -ratelaws(425) -ratelaws(426);
    Dspecies(64) = ratelaws(193) -ratelaws(207) +ratelaws(208) -ratelaws(211) +ratelaws(212);
    Dspecies(65) = ratelaws(194) -ratelaws(205) +ratelaws(206) -ratelaws(209) +ratelaws(210);
    Dspecies(66) = ratelaws(195) +ratelaws(205) -ratelaws(206) -ratelaws(213) +ratelaws(318);
    Dspecies(67) = ratelaws(196) +ratelaws(207) -ratelaws(208) -ratelaws(214) +ratelaws(319);
    Dspecies(68) = ratelaws(197) +ratelaws(209) -ratelaws(210) -ratelaws(223) +ratelaws(322);
    Dspecies(69) = ratelaws(198) +ratelaws(211) -ratelaws(212) -ratelaws(224) +ratelaws(323);
    Dspecies(70) = ratelaws(199) +ratelaws(200) +ratelaws(201) +ratelaws(202) +ratelaws(203) +ratelaws(204) +ratelaws(313) +ratelaws(314) +ratelaws(315) +ratelaws(316) -ratelaws(317) +ratelaws(396) +ratelaws(397) +ratelaws(398) +ratelaws(433) +ratelaws(434);
    Dspecies(71) = ratelaws(213) +ratelaws(309) -ratelaws(318) -ratelaws(320);
    Dspecies(72) = ratelaws(214) +ratelaws(310) -ratelaws(319) -ratelaws(321);
    Dspecies(73) = ratelaws(219) +ratelaws(220);
    Dspecies(74) = ratelaws(221) +ratelaws(269) +ratelaws(270) +ratelaws(271) +ratelaws(272) +ratelaws(273) +ratelaws(274) +ratelaws(275) +ratelaws(276) +ratelaws(277) +ratelaws(278) +ratelaws(354) +ratelaws(355) +ratelaws(356) -ratelaws(392) +ratelaws(419) +ratelaws(420);
    Dspecies(75) = ratelaws(222) +ratelaws(279) +ratelaws(280) +ratelaws(281) +ratelaws(282) +ratelaws(283) +ratelaws(284) +ratelaws(285) +ratelaws(286) +ratelaws(287) +ratelaws(288) +ratelaws(357) +ratelaws(358) +ratelaws(359) -ratelaws(393) +ratelaws(421) +ratelaws(422);
    Dspecies(76) = ratelaws(223) +ratelaws(311) -ratelaws(322) -ratelaws(366) -ratelaws(367) -ratelaws(368) -ratelaws(369) -ratelaws(370) -ratelaws(371) -ratelaws(372) -ratelaws(373) -ratelaws(374) -ratelaws(375) -ratelaws(376) -ratelaws(377) -ratelaws(378) -ratelaws(427) -ratelaws(428);
    Dspecies(77) = ratelaws(224) +ratelaws(312) -ratelaws(323) -ratelaws(379) -ratelaws(380) -ratelaws(381) -ratelaws(382) -ratelaws(383) -ratelaws(384) -ratelaws(385) -ratelaws(386) -ratelaws(387) -ratelaws(388) -ratelaws(389) -ratelaws(390) -ratelaws(391) -ratelaws(429) -ratelaws(430);
    Dspecies(78) = ratelaws(289) +ratelaws(290) +ratelaws(291) +ratelaws(292) +ratelaws(293) +ratelaws(294) +ratelaws(295) +ratelaws(296) +ratelaws(297) +ratelaws(298) +ratelaws(299) +ratelaws(300) +ratelaws(301) +ratelaws(302) +ratelaws(303) +ratelaws(304) +ratelaws(305) +ratelaws(306) +ratelaws(307) +ratelaws(308) +ratelaws(360) +ratelaws(361) +ratelaws(362) +ratelaws(363) +ratelaws(364) +ratelaws(365) -ratelaws(394) -ratelaws(395) +ratelaws(423) +ratelaws(424) +ratelaws(425) +ratelaws(426);
    Dspecies(79) = ratelaws(320) +ratelaws(392);
    Dspecies(80) = ratelaws(321) +ratelaws(393);
    Dspecies(81) = ratelaws(366) +ratelaws(367) +ratelaws(368) +ratelaws(369) +ratelaws(370) +ratelaws(371) +ratelaws(372) +ratelaws(373) +ratelaws(374) +ratelaws(375) +ratelaws(376) +ratelaws(377) +ratelaws(378) +ratelaws(395) +ratelaws(427) +ratelaws(428) -ratelaws(431);
    Dspecies(82) = ratelaws(379) +ratelaws(380) +ratelaws(381) +ratelaws(382) +ratelaws(383) +ratelaws(384) +ratelaws(385) +ratelaws(386) +ratelaws(387) +ratelaws(388) +ratelaws(389) +ratelaws(390) +ratelaws(391) +ratelaws(394) +ratelaws(429) +ratelaws(430) -ratelaws(432);
    Dspecies(83) = ratelaws(431) +ratelaws(432);

end


end
