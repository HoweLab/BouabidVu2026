
% combine individual mouse fib loc tables into a cohort table
% Mai-Anh Vu, 2024/09/23
function fib = cohort_fib_table(data_dir,mice)
    for m = 1:numel(mice)
        mouse = mice{m};
        tbl = readtable(fullfile(data_dir,mouse,'fiber_table.xlsx'));
        tbl = tbl(tbl.included,:);
        tbl.mouse = repmat({mouse},size(tbl,1),1);
        tbl.ROI_orig = tbl.ROI;
        % concatenate
        if m == 1 % if it's the first mouse, start the table
            fib = tbl;        
        else % for subsequent mice, concatenate 
            fib = cat(1,fib,tbl);
        end    
    end
    fib.ROI = vec(1:size(fib,1));
end
    
