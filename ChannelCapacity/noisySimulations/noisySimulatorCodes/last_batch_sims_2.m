ant_val = [0.0, 0.507471, 4.1422, 6.8024, 9.5051, 12.2511, 15.2534, 18.9914,...
            70.7456];
sigma   = [0.25, 0.5, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0];

for idx=1:14
    for jdx=1:9
        noisySimulator_extend(ant_val(idx), sigma(jdx));
    end
end