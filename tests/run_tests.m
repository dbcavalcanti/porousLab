function run_tests()
    clear; close all; clc;

    test_dir = fullfile(pwd, 'files');
    ref_dir  = fullfile(pwd, 'references');

    test_files = dir(fullfile(test_dir, 'test_*.m'));
    tol = 1e-8;

    total = 0;
    passed = 0;
    failed = 0;
    skipped = 0;

    for k = 1:length(test_files)
        test_name = test_files(k).name(1:end-2);  % strip .m
        fprintf('Running %s...\n', test_name);

        clear mdl;
        run(fullfile(test_dir, test_files(k).name));
        total = total + 1;

        ref_file = fullfile(ref_dir, [test_name '_ref.mat']);
        if ~exist(ref_file, 'file')
            warning('⚠️ Reference file missing: %s\n', ref_file);
            skipped = skipped + 1;
            continue;
        end

        load(ref_file, 'U_ref');

        if ~isequal(size(mdl.U), size(U_ref))
            fprintf('❌ %s failed: size mismatch (mdl.U is [%s], U_ref is [%s])\n\n', ...
                test_name, num2str(size(mdl.U)), num2str(size(U_ref)));
            failed = failed + 1;
            continue;
        end

        err = norm(mdl.U - U_ref) / norm(U_ref);

        if err < tol
            fprintf('✅ %s passed (rel. error = %.2e)\n\n', test_name, err);
            passed = passed + 1;
        else
            fprintf('❌ %s failed (rel. error = %.2e)\n\n', test_name, err);
            failed = failed + 1;
        end
    end

    % Final summary
    fprintf('======== Test Summary ========\n');
    fprintf('Total tests:     %d\n', total);
    fprintf('Passed:          %d\n', passed);
    fprintf('Failed:          %d\n', failed);
    fprintf('Skipped (no ref):%d\n', skipped);
    fprintf('===============================\n');
end
