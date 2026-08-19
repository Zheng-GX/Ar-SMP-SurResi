function samples = TumbleAngleSampling1(n_samples)
    % 
    
    % 
    bins = [-180, -160, -126, -92, -58, -24, 10, 44, 78, 112, 146, 180];
    densities = [0, 0.00155972, 5.03e-04, 0.00208786, 9.87e-04, 0.00728036, ...
                 0.00939783, 0.00310464, 0.00209279, 0.00100691, 0.00157947];
    
    max_density = max(densities);
    samples = zeros(n_samples, 1);
    count = 0;
    
    while count < n_samples
        % 
        candidate = -180 + 360 * rand();
        
        % 
        for i = 1:length(bins)-1
            if candidate >= bins(i) && candidate < bins(i+1)
                density_val = densities(i);
                break;
            end
        end
        
        % 
        if rand() <= density_val / max_density
            count = count + 1;
            samples(count) = candidate;
        end
    end
    samples=samples*pi/180;
end