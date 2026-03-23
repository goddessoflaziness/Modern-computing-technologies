import Pkg
Pkg.activate(dirname(@__DIR__))

using Dirichlet_problem_solver
using LinearAlgebra

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

    print("Введите количество шагов разбиения сетки: ")
    n_step = parse(Int, readline())

    x, u_num = Solve_Laplace_eq(f, n_step, a, b)

    u_ex = u_exact.(x)

    err_C = maximum(abs.(u_num - u_ex))
    err_L2 = norm(u_num - u_ex)
    println("Количество шагов = $n_step, шаг разбиения = $(1/n_step)")
    println("Ошибка в C-норме = $err_C")
    println("Ошибка в L2-норме = $err_L2")
end

main()