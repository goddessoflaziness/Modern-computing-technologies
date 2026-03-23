"""
Tridiagonal_matrix_algorithm(
                a::Vector, b::Vector, c::Vector, f::Vector)

Решает систему линейных уравнений с трехдиагональной матрицей
размера n*n методом прогонки

# Аргументы
- 'a' : вектор длины n-1, содержащий поддиагональные элементы матрицы,
- 'b' : вектор длины n, содержащий диагональные элементы матрицы,
- 'c' : вектор длины n-1, содержащий наддиагональные элементы матрицы,
- 'f' : вектор - правая часть системы.

# Возвращаемое значение
u – массив, содержащий решение системы

```jldoctest
julia> a = [-1.0, -1.0]; b = [2.0, 2.0, 2.0]; c = [-1.0, -1.0]; f = [1.0, 1.0, 1.0];
julia> Tridiagonal_matrix_algorithm(a, b, c, f)
3-element Vector{Float64}:
 1.5
 2.0
 1.5
 ```
 """
function Tridiagonal_matrix_algorithm(a::Vector, b::Vector, c::Vector, f::Vector)
    n = length(b)
    coef = similar(b, n)
    u = zeros(n)

    coef[1] = c[1] / b[1]
    u[1] = f[1] / b[1]

    for idx = 2:n
        denominator = b[idx] - a[idx - 1] * coef[idx - 1]
        if (idx < n) 
            coef[idx] = c[idx] / denominator
        end
        u[idx] = (f[idx] - a[idx - 1] * u[idx - 1]) / denominator
    end

    for idx = n - 1:-1:1
        u[idx] -= coef[idx] * u[idx + 1]
    end
    return u
end

function Generate_rhs(func::Function, n_steps)
    h = 1 / n_steps
    x = [i * h for i = 0:(n_steps)]
    h = h^2
    return x, [func(x[i]) * h for i = 1:(n_steps - 1)]
end

"""
Solve_Laplace_eq(func::Function, n_step, a=0.0, b=0.0)

Численно решает краевую задачу Дирихне 
для уравнения Лапласа на интревале (0; 1)

# Аргументы
- 'func' : функция - правая часть уравнения Лапласа
- 'n_step' : количество шагов разбиения сетки (всего n_step + 1 узлов)
- 'a' : значение на левой границе
- 'b' : значение на правой границе

# Возвращаемое значение
Кортеж, где
x - значение узлов сетки
res – массив, содержащий решение в узлах сетки
"""
function Solve_Laplace_eq(func::Function, n_step, a=0.0, b=0.0)
    a_diag = fill(-1.0, n_step - 2)
    b_diag = fill(2.0, n_step - 1)
    c_diag = fill(-1.0, n_step - 2)
    
    x, f = Generate_rhs(func, n_step)
    f[1] += a
    f[end] += b
    
    res = zeros(n_step + 1)
    res[2:end - 1] = Tridiagonal_matrix_algorithm(a_diag, b_diag, c_diag, f)
    res[1] = a
    res[end] = b
    return x, res
end