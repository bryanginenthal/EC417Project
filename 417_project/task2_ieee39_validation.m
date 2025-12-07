function task2_ieee39_validation()
    fprintf('IEEE-39 validation\n')

    mpc = case39;

    nb = size(mpc.bus,1);
    ng = size(mpc.gen,1);
    nl = size(mpc.branch,1);

    fprintf('buses: %d\n', nb);
    fprintf('generators: %d\n', ng);
    fprintf('branches: %d\n', nl);

    bt = mpc.bus(:,2);
    nPQ  = sum(bt == 1);
    nPV  = sum(bt == 2);
    nREF = sum(bt == 3);

    fprintf('PQ buses: %d\n', nPQ);
    fprintf('PV buses: %d\n', nPV);
    fprintf('slack buses: %d\n', nREF);

    Pg  = mpc.gen(:,2);
    Qg  = mpc.gen(:,3);
    Qmax = mpc.gen(:,4);
    Qmin = mpc.gen(:,5);
    Pmax = mpc.gen(:,9);
    Pmin = mpc.gen(:,10);
    gb   = mpc.gen(:,1);

    fprintf('generator limit issues (if any):\n')
    for k = 1:ng
        p_hi = Pg(k) > Pmax(k) + 1e-6;
        p_lo = Pg(k) < Pmin(k) - 1e-6;
        q_hi = Qg(k) > Qmax(k) + 1e-6;
        q_lo = Qg(k) < Qmin(k) - 1e-6;
        if p_hi || p_lo || q_hi || q_lo
            fprintf('  gen @ bus %2d: Pg=%.2f [%.2f %.2f], Qg=%.2f [%.2f %.2f]\n', ...
                gb(k), Pg(k), Pmin(k), Pmax(k), Qg(k), Qmin(k), Qmax(k));
        end
    end

    [Ybus, ~, ~] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
    Yfull = full(Ybus);
    bad = any(isnan(Yfull(:))) || any(isinf(Yfull(:)));
    if bad
        fprintf('Ybus has NaN/Inf entries\n')
    else
        fprintf('Ybus ok (%dx%d)\n', size(Ybus,1), size(Ybus,2));
    end

    Ym = abs(Ybus);
    Ym(1:size(Ym,1)+1:end) = 0;
    coupled = sum(Ym,2);
    iso = find(coupled < 1e-8);
    if isempty(iso)
        fprintf('no isolated buses\n')
    else
        fprintf('isolated buses: %s\n', mat2str(iso'))
    end

    if exist('mpoption','file')
        mpopt = mpoption('pf.alg','NR','verbose',1,'out.all',0);
        [res, ok] = runpf(mpc, mpopt);
    else
        [res, ok] = runpf(mpc);
    end

    if ~ok
        fprintf('AC power flow did not converge\n')
        return
    end

    fprintf('AC power flow converged\n')

    Vm = res.bus(:,8);
    Va = res.bus(:,9);

    fprintf('Vm min/max: [%.4f %.4f]\n', min(Vm), max(Vm));
    fprintf('Va min/max: [%.2f %.2f] deg\n', min(Va), max(Va));

    fprintf('first few generator outputs:\n')
    for k = 1:min(5, size(res.gen,1))
        fprintf('  gen @ bus %2d: Pg=%.2f, Qg=%.2f\n', ...
            res.gen(k,1), res.gen(k,2), res.gen(k,3));
    end
end
