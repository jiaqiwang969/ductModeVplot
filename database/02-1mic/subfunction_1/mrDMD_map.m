function [map, low_f_cutoff] = mrDMD_map(mrdmd)
% function [map, low_f_cutoff] = mrDMD_map(mrdmd)

[levels, M] = size(mrdmd);

map = zeros(levels, M);

low_f_cutoff = zeros(levels+1, 1);
for i = 1:levels,
    chunks = 2^(i-1);
    K = M / chunks;
    for j = 1:chunks,
        f = abs(imag(mrdmd{i, j}.omega));
        P = mrdmd{i, j}.P;

        P = P(f >= low_f_cutoff(i));
        if ~isempty(P),
            map(i, (1:K)+(j-1)*K) = mean(P);
        end;
    end;
    
    low_f_cutoff(i+1) = mrdmd{i,1}.rho;
end;