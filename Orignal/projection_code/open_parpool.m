function open_parpool(no_par_processes)

% parallel processes? set this to zero for just one process
if nargin<1
    no_par_processes=16; 
end

disp(['PPN: ',num2str(no_par_processes)]);

cp=gcp('nocreate');
if ~isempty(cp)
    if  cp.NumWorkers~=no_par_processes
        delete(cp);
        if no_par_processes>0
            parpool(no_par_processes);
        end
    end
else
    if no_par_processes>0
        parpool(no_par_processes);
    end
end
