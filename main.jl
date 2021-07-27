include("hamiltonians.jl")
include("observables.jl")

using DifferentialEquations, Plots, JLD

function main_approx(n_cut, r_end, points)
    # set parameters:
    p = (0.1, n_cut);  # p = (r, n_cut). r will be updated in the loop below
    timespan = (0.0, 8.0);  # integration time span
    r_vals = range(0.001, r_end, length=points);  # range and number of r values to be calculated

    # initial conditions for the approximate Hamiltonian: 
    u0_approx = zeros(ComplexF64, p[2]+2, p[2]+2, p[2]+2, p[2]+2);
    # initializes an array with one additional dimension to avoid out-of-bounds indexing error
    # u[n_cut+2, :, :, :] etc. will all be zero

    u0_approx[1, 1, 1, 1] = 1 + 0im;  # start with 0 atoms in all modes

    # Preallocating vectors for the results:
    approach1_vals = zeros(length(r_vals));
    approach2_vals = zeros(length(r_vals));
    approach2_gamma09_vals = zeros(length(r_vals));
    approach2_gamma07_vals = zeros(length(r_vals));
    approach3_vals = zeros(length(r_vals));
    approach3_gamma09_vals = zeros(length(r_vals));
    approach3_gamma07_vals = zeros(length(r_vals));

    loss = zeros(n_cut+1, n_cut+1);

    #start = time(); # optional tracking of the execution time (uncomment the corresponding lines in the loop)
    index = 0;

    # Solving the system and calculating observables:
    for r in r_vals
        p = (r, n_cut);
        prob_approx = ODEProblem(H_approx!, u0_approx, timespan, p);  
        if r < 0.05  # for very small r, the tolerance needs to be smaller to avoid numerical errors
            sol = solve(prob_approx, Tsit5(), abstol=1e-9, reltol=1e-9, saveat=[5.0,6.0,7.0,8.0]);
        else
            sol = solve(prob_approx, Tsit5(), abstol=1e-8, reltol=1e-8, saveat=[5.0,6.0,7.0,8.0]);
        end
        
        index += 1;
        
        approach1_vals[index] = approach1(sol, p[2]);
        app2, app3 = approach23(sol, p[2]);
        app2_09, app3_09 = approach23_loss!(sol, p[2], 0.9, loss);
        app2_07, app3_07 = approach23_loss!(sol, p[2], 0.7, loss);
        approach2_vals[index] = app2;
        approach3_vals[index] = app3;
        approach2_gamma09_vals[index] = app2_09;
        approach3_gamma09_vals[index] = app3_09;
        approach2_gamma07_vals[index] = app2_07;
        approach3_gamma07_vals[index] = app3_07;
    end

    #elapsed = time() - start;
    #println("elapsed time:", elapsed);

    # save the data for later use:
    save("data_approx.jld", "approach1_vals", approach1_vals,
        "approach2_vals", approach2_vals,
        "approach3_vals", approach3_vals,
        "approach3_gamma07_vals", approach3_gamma07_vals,
        "approach3_gamma09_vals", approach3_gamma09_vals,
        "approach2_gamma07_vals", approach2_gamma07_vals,
        "approach2_gamma09_vals", approach2_gamma09_vals)

    # Plotting:
    p1_approx = plot(r_vals, approach1_vals, label="(I)");
    plot!(p1_approx, r_vals, approach3_vals, label="(III) γ=1.0");
    plot!(p1_approx, r_vals, approach3_gamma09_vals, label="(III) γ=0.9");
    plot!(p1_approx, r_vals, approach3_gamma07_vals, label="(III) γ=0.7");
    plot!(xlims = (0, r_end), ylims = (1.95, 2.9), xticks = 0:0.1:r_end, yticks = 2.0:0.1:2.8, xlabel = "r", ylabel="B");
    plot!(xtickfontsize=16, ytickfontsize=16, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, legend=:bottomleft)
    savefig(p1_approx, "p1_approx.pdf")

    p2_approx = plot(r_vals, approach2_vals, label="(II) γ=1.0");
    plot!(p2_approx, r_vals, approach2_gamma09_vals, label="(II) γ=0.9");
    plot!(p2_approx, r_vals, approach2_gamma07_vals, label="(II) γ=0.7");
    plot!(xlims = (0, r_end), ylims = (1.75, 2.35), xticks = 0:0.1:r_end, yticks = 1.8:0.1:2.3, xlabel = "r", ylabel="B");
    plot!(xtickfontsize=16, ytickfontsize=16, legendfontsize=16, xguidefontsize=16, yguidefontsize=16, legend=:topleft)
    savefig(p2_approx, "p2_approx.pdf")
end


function main_exact(N, n_cut, r_end, points)

    # set parameters:
    p = (0.1, n_cut, N); # p = (r, n_cut, N). r will be updated in the loop below
    timespan = (0.0, 8.0);  # integration time span
    r_vals = range(0.003, r_end, length=points);  # range and number of r values to be calculated

    # initial conditions for the exact  Hamiltonian:
    u0_exact = zeros(ComplexF64, p[2]+2, p[2]+2, p[3]+3, p[2]+2, p[2]+2, p[3]+3);
    u0_exact[1, 1, p[3]+1, 1, 1, 1] = 1 + 0im;  # start with N atoms in k3 mode (zero-mode of site A)

    # Preallocating vectors for the results:
    approach1_vals = zeros(length(r_vals));
    approach2_vals = zeros(length(r_vals));
    approach2_gamma09_vals = zeros(length(r_vals));
    approach2_gamma07_vals = zeros(length(r_vals));
    approach2_gamma095_vals = zeros(length(r_vals));
    approach3_vals = zeros(length(r_vals));
    approach3_gamma09_vals = zeros(length(r_vals));
    approach3_gamma07_vals = zeros(length(r_vals));

    loss = zeros(n_cut+1, n_cut+1);

    #start = time(); # optional tracking the execution time (uncomment the corresponding lines in the loop)
    index = 0;

    for r in r_vals
        p = (r, n_cut, N);
        prob_exact = ODEProblem(H_exact!, u0_exact, timespan, p);

        #checkpoint = time()  # optional time tracking
        if r < 0.05  # for very small r, the tolerance needs to be smaller to avoid numerical errors
            sol = solve(prob_exact, Tsit5(), abstol=1e-9, reltol=1e-9, saveat=[5.0,6.0,7.0,8.0]);
        else
            sol = solve(prob_exact, Tsit5(), abstol=1e-8, reltol=1e-8, saveat=[5.0,6.0,7.0,8.0]);
        end
        
        index += 1;
        
        println("Solved ", index)
        #elapsed = time() - checkpoint;
        #println("Time elapsed for ", index, " :", elapsed);
        
        approach1_vals[index] = approach1_exact(sol, p[2], p[3]);
        app2, app3 = approach23_exact(sol, p[2], p[3]);
        app2_09, app3_09 = approach23_loss_exact!(sol, p[2], p[3], 0.9, loss);
        app2_07, app3_07 = approach23_loss_exact!(sol, p[2], p[3], 0.7, loss);
        app2_095, app3_095 = approach23_loss_exact!(sol, p[2], p[3], 0.95, loss);
        approach2_vals[index] = app2;
        approach3_vals[index] = app3;
        approach2_gamma09_vals[index] = app2_09;
        approach3_gamma09_vals[index] = app3_09;
        approach2_gamma07_vals[index] = app2_07;
        approach3_gamma07_vals[index] = app3_07;
        approach2_gamma095_vals[index] = app2_095;
        println("Done with ", index)
    end

    #elapsed = time() - start;
    #println("elapsed time:", elapsed);

    # save the data for later use:
    save("data_exact.jld", "approach1_vals", approach1_vals,
        "approach2_vals", approach2_vals,
        "approach3_vals", approach3_vals,
        "approach3_gamma07_vals", approach3_gamma07_vals,
        "approach3_gamma09_vals", approach3_gamma09_vals,
        "approach2_gamma07_vals", approach2_gamma07_vals,
        "approach2_gamma09_vals", approach2_gamma09_vals,
        "approach2_gamma095_vals", approach2_gamma095_vals)

    # Plotting:
    p1 = plot(r_vals, approach1_vals, label="(I)");
    plot!(p1, r_vals, approach3_vals, label="(III) γ=1.0");
    plot!(p1, r_vals, approach3_gamma09_vals, label="(III) γ=0.9");
    plot!(p1, r_vals, approach3_gamma07_vals, label="(III) γ=0.7");
    plot!(xlims = (0, r_end), ylims = (1.95, 2.9), xticks = 0:0.1:r_end, yticks = 2.0:0.1:2.8, xlabel = "r", ylabel="B");
    plot!(xtickfontsize=14, ytickfontsize=14, legendfontsize=14, xguidefontsize=14, yguidefontsize=14, legend=:bottomleft)
    savefig(p1, "p1_exact.pdf")

    p2 = plot(r_vals, approach2_vals, label="(II) γ=1.0");
    plot!(p2, r_vals, approach2_gamma095_vals, label="(II) γ=0.95");
    plot!(p2, r_vals, approach2_gamma09_vals, label="(II) γ=0.9");
    plot!(p2, r_vals, approach2_gamma07_vals, label="(II) γ=0.7");
    plot!(xlims = (0, r_end), ylims = (1.75, 2.35), xticks = 0:0.1:r_end, yticks = 1.8:0.1:2.3, xlabel = "r", ylabel="B");
    plot!(xtickfontsize=14, ytickfontsize=14, legendfontsize=14, xguidefontsize=14, yguidefontsize=14, legend=:topleft)
    savefig(p2, "p2_exact.pdf")
end

# Examples for the execution of the code (uncomment and set parameters to desired values):
#   main_approx(25, 1.0, 15)
#   main_exact(20, 10, 0.5, 10)
