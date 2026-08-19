function samples = TumbleAngleSampling_Berg(n_samples)
    % 
    
    % 
    % AngleData=[5 	15 	25 	35 	45 	55 	65 	75 	85 	95 	105 	115 	125 	135 	145 	155 	165 	175];
bins = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160  170 180];
densities = [20 	54 	72 	121 	142 	139 	135 	93 	90 	57 	52 	54 	41 	41 	18 	12 	10 	4 ];

    
    max_density = max(densities);
    samples = zeros(n_samples, 1);
    count = 0;
    
    while count < n_samples
        % 
        candidate = 180* rand();
        
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
            samples(count) = candidate*pi/180;
        end
    end
end