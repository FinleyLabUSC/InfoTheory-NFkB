ant_val = [0.0, 0.7, 4.0, 4.5, 5.5, 7.5, 9.0, 11.0,...
            13.0, 14.5, 16.5, 18.0, 19.0, 70.0];
sigma   = [0.25, 0.5, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 5.0];

for idx=1:14
    for jdx=1:9
        noisySimulator_extend(ant_val(idx), sigma(jdx));
    end
end