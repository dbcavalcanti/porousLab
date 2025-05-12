function run_tests()
    clear; close all; clc;

    test_dir = fullfile(pwd, 'files');
    ref_dir  = fullfile(pwd, 'references');
    
    test_files = dir(fullfile(test_dir, 'test_*.m'));
    tol = 1e-8;

    for k = 1:length(test_files)
        test_name = test_files(k).name(1:end-2);  % strip .m
        fprintf('Running %s...\n', test_name);

        clear mdl;
        run(fullfile(test_dir, test_files(k).name));

        ref_file = fullfile(ref_dir, [test_name '_ref.mat']);
        if ~exist(ref_file, 'file')
            warning('Reference file missing: %s', ref_file);
            continue;
        end

        load(ref_file, 'U_ref');
        err = norm(mdl.U - U_ref) / norm(U_ref);

        if err < tol
            fprintf('✅ %s passed (rel. error = %.2e)\n\n', test_name, err);
        else
            fprintf('❌ %s failed (rel. error = %.2e)\n\n', test_name, err);
        end
    end
end