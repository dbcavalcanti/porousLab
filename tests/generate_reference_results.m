function generate_reference_results()
    clear; close all; clc;

    % Base directory
    base_dir = fileparts(mfilename('fullpath'));

    % Interest paths
    src_dir = fullfile(base_dir,  '..', 'src');
    files_dir  = fullfile(base_dir,'files');
    ref_dir  = fullfile(base_dir,'references');

    % Add path to the source code
    addpath(genpath(src_dir));

    % Check existence of the folders
    if ~exist(files_dir, 'dir')
        error('Missing files folder');
    end
    if ~exist(ref_dir, 'dir')
        mkdir(ref_dir);
    end

    % Search for all test_*.m scripts recursively
    files = dir(fullfile(base_dir, '**', 'test_*.m'));

    for i = 1:length(files)
        file = files(i);
        [~, test_name] = fileparts(file.name);
        script_path = fullfile(file.folder, file.name);
        fprintf('Running %s...\n', script_path);

        try
            clear mdl U_ref;
            run(script_path);

            if ~exist('mdl', 'var') || ~isprop(mdl, 'U') || isempty(mdl.U)
                warning('Skipped %s: mdl.U not found.', test_name);
                continue;
            end

            % Save reference result
            U_ref = mdl.U;
            save(fullfile(ref_dir, [test_name '_ref.mat']), 'U_ref');
            fprintf('Saved reference: %s\n\n', [test_name '_ref.mat']);

        catch ME
            warning('Failed to run %s: %s', test_name, ME.message);
        end
    end
end
