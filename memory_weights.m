function weights = memory_weights(sigma_record, sigma_gen, memory_size, memory_sf)
%MEMORY_WEIGHTS 
%   此处显示详细说明
    
    
    if numel(sigma_record) < sigma_gen * 2 || numel(sigma_record) == 0
        weights = ones(memory_size,1) / memory_size;
    else
        max_pos = 0.3;
        last_sigma_mean= mean(sigma_record(end - sigma_gen*2+1:end-sigma_gen));
        cur_sigma_mean = mean(sigma_record(end - sigma_gen+1:end));
        weights = ones(memory_size,1);
        if last_sigma_mean < cur_sigma_mean
            % sigma increases, scale factor F increases
            [~,index] = max(memory_sf);
            weights(1:end) = (1-max_pos)/(memory_size-1);
            weights(index) = 0.5;
        else
            % sigma decreases, scale factor F decreases
            [~,index] = min(memory_sf);
            weights(1:end) = (1-max_pos)/(memory_size-1);
            weights(index) = 0.5;
        end
    end
    
    for i = 2 : memory_size
        weights(i) = weights(i) + weights(i-1);
    end
    
end

