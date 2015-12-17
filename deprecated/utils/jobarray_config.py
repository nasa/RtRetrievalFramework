class config:
    tc_dir = '/home/mcduffie/oco_l2/auto_testing/testcases'
    output_dir = '/home/mcduffie/oco_l2/auto_testing/run_output'

    max_run_time = 180
    queue_name = 'short'
    binary_name = '/scratch05/mcduffie/L2_EXE/bin/oco_l2.f90'
    num_simul_jobs = 255

    local_dir = '/lscratch/mcduffie'

    cp_bin = '/bin/cp'
    mv_bin = '/bin/mv'   
    rm_bin = '/bin/rm'

    run_dir_search_glob = 'run.*'
    testcase_dir_name_re = 'testcase_[^/]+'
