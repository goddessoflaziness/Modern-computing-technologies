import Pkg
Pkg.activate(dirname(@__DIR__))

using CSV, DataFrames, Plots

# Чтение данных из result.txt
function parse_result_txt(filename)
    sizes = Int[]
    iterations = Int[]
    prec_times = Float64[]
    iter_times = Float64[]
    l2_errors = Float64[]
    c_errors = Float64[]
    
    open(filename, "r") do f
        lines = readlines(f)
        for line in lines
            if occursin("Size of matrix:", line)
                push!(sizes, parse(Int, split(line, ':')[2]))
            elseif occursin("Number of iterations:", line)
                push!(iterations, parse(Int, split(line, ':')[2]))
            elseif occursin("Preconditioner Time:", line)
                push!(prec_times, parse(Float64, split(line, ':')[2]))
            elseif occursin("Iterations time:", line)
                push!(iter_times, parse(Float64, split(line, ':')[2]))
            elseif occursin("L2-norm error:", line)
                push!(l2_errors, parse(Float64, split(line, ':')[2]))
            elseif occursin("C-norm error:", line)
                push!(c_errors, parse(Float64, split(line, ':')[2]))
            end
        end
    end

    return DataFrame(size = sizes, iterations = iterations,
                     prec_time = prec_times, iter_time = iter_times,
                     l2_error = l2_errors, c_error = c_errors)
end

# Чтение данных из results.csv
function parse_results_csv(filename)
    df = CSV.read(filename, DataFrame; header = false)
    rename!(df, [:size, :h, :c_error, :l2_error, :prec_time, :iter_time])
    return df
end

function main()
    df_txt = parse_result_txt("result.txt")
    df_csv = parse_results_csv("results.csv")

    df_csv[!, :iterations] = df_txt.iterations

    df_csv[!, :total_time] = df_csv.prec_time + df_csv.iter_time

    # График стравнения ошибок
    plot1 = plot(df_csv.iterations, df_csv.l2_error,
                xscale = :log10, yscale = :log10,
                label = "L2‑норма ошибки",
                marker = :circle, linewidth = 2,
                xlabel = "Число итераций", ylabel = "Ошибка",
                title = "Зависимость ошибки от числа итераций",
                legend = :topright)
    plot!(df_csv.iterations, df_csv.c_error,
        label = "C‑норма ошибки",
        marker = :circle, linewidth = 2)

    # График сравнения времени вычисления
    plot2 = plot(df_csv.iterations, df_csv.total_time,
                xscale = :log10, yscale = :log10,
                label = "Общее время (предобусл. + итерации)",
                marker = :circle, linewidth = 2,
                xlabel = "Число итераций", ylabel = "Время (с)",
                title = "Зависимость времени вычислений от числа итераций",
                legend = :topleft)
    plot!(df_csv.iterations, df_csv.prec_time,
        label = "Время предобуславливателя",
        marker = :circle, linewidth = 2)
    plot!(df_csv.iterations, df_csv.iter_time,
    label = "Время итераций",
    marker = :circle, linewidth = 2)

    savefig(plot1, "errors.png")
    savefig(plot2, "time.png")
end

main()
