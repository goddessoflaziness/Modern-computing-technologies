import Pkg
Pkg.activate(dirname(@__DIR__))

using Dirichlet_problem_solver
using LinearAlgebra
using Plots

function u_exact(x)
    return sin(3x) * cos(4x)
end

# -u'' получено аналитически
function f(x)
    return 25 * sin(3x) * cos(4x) + 24 * cos(3x) * sin(4x)
end

function main()
    a = u_exact(0.0)
    b = u_exact(1.0)

    n_steps_list = [10, 16, 20, 50, 100, 200, 500, 1000, 5000, 10000, 100000]
    errors_С = []
    errors_L2 = []

    for n_step in n_steps_list

        x, u_num = Solve_Laplace_eq(f, n_step, a, b)

        u_ex = u_exact.(x)

        err_C = maximum(abs.(u_num - u_ex))
        err_L2 = sqrt(sum((1.0 / n_step) * (abs.(u_num - u_ex)).^2))
        push!(errors_С, err_C)
        push!(errors_L2, err_L2)

        if (n_step == 16)
            p1 = plot(x, [u_num, u_ex],
                label=["Численное решение" "Аналитическое решение"],
                xlabel="x", ylabel="u(x)",
                title="Сравнение решений, n = $n_step",
                linestyle=[:solid :solid], linewidth=2,
                legend=:topleft)

            l = @layout [a; b]
            plot(p1, layout=l, size=(800, 1000))
            savefig("report_plots_$n_step.png")
        end
    end

    h = [1/n for n in n_steps_list]
    h² = [(1/n)^2 for n in n_steps_list]

    p2 = plot(n_steps_list, [errors_С, errors_L2, h, h²],
        label=["Ошибка в C-норме" "Ошибка в L2-норме" "h" "h²"],
        xscale=:log10, yscale=:log10,
        xlabel="Количество шагов", ylabel="Ошибка",
        title="Зависимость ошибок от числа узлов сетки",
        linestyle=[:solid, :solid], linewidth=4,
        legend=:topleft)

    l = @layout [a; b]
    plot(p2, layout=l, size=(800, 1000))
    savefig("report_errors.png")
end

main()