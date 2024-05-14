// Create the selected ODE system solver
dictionary odeDict;
odeDict.add("solver", args.getOrDefault<word>("ODESolver", "RKDP45"));
motionODE ode(6.0);
autoPtr<ODESolver> motionSolver(ODESolver::New(ode, odeDict));
motionSolver->relTol() = 1e-9;

scalar tStart(0.0);
scalarField sStart(ode.nEqns());