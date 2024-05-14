// Disable dynamic teams explicitly
omp_set_dynamic(0);

if (args.found("np"))
{
    label threadNum(args.get<label>("np"));

    omp_set_num_threads
    (
        (threadNum < omp_get_num_procs())
        ?
        threadNum
        :
        omp_get_num_procs() - 1
    );

    Info<< nl
        << "Running with "
        << omp_get_max_threads() << "/" << omp_get_num_procs()
        << " threads"
        << nl << endl;
}
else
{
    omp_set_num_threads(0);

    Info << nl << "Running in serial mode" << nl << endl;
}
