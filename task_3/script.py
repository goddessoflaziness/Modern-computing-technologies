import csv

"""
GEO-скрипт для Gmsh на основе координат из CSV.
"""
def generate_geo_script(csv_filename, output_filename, refine_point_x, refine_point_y, refine_size, boundary_size):
    points = []
    try:
        with open(csv_filename, 'r', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)
            for i, row in enumerate(reader):
                x, y = float(row[0]), float(row[1])
                points.append((x, y))
    except FileNotFoundError:
        print(f"Ошибка: Файл '{csv_filename}' не найден.")
        return
    except Exception as e:
        print(f"Ошибка при чтении CSV: {e}")
        return

    if not points:
        print("Ошибка: Не удалось загрузить точки из CSV.")
        return

    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write('// =====================================================\n')
        f.write('// Сетка для Московской области\n')
        f.write('// =====================================================\n\n')

        f.write('SetFactory("OpenCASCADE");\n\n')

        f.write('// Точки границы\n')
        for idx, (x, y) in enumerate(points):
            tag = idx + 1
            f.write(f'Point({tag}) = {{{x}, {y}, 0, {boundary_size}}};\n')
        f.write('\n')

        f.write('// Линии границы\n')
        num_points = len(points)
        for i in range(num_points):
            start_tag = i + 1
            end_tag = (i + 1) % num_points + 1
            line_tag = i + 1
            f.write(f'Line({line_tag}) = {{{start_tag}, {end_tag}}};\n')
        f.write('\n')

        f.write('// Контур и поверхность\n')
        line_loop_tags = ', '.join(str(i+1) for i in range(num_points))
        f.write(f'Curve Loop(1) = {{{line_loop_tags}}};\n')
        f.write(f'Plane Surface(1) = {{1}};\n\n')

        f.write('// Точка сгущения\n')
        f.write(f'Point({num_points+1}) = {{{refine_point_x}, {refine_point_y}, 0, {refine_size}}};\n')
        f.write(f'Point{{{num_points+1}}} In Surface{{1}};\n\n')

        f.write('// Настройка полей для плавного перехода размера ячеек\n')
        f.write('Field[1] = Distance;\n')
        f.write(f'Field[1].PointsList = {{{num_points+1}}};\n\n')
        f.write('Field[2] = Threshold;\n')
        f.write('Field[2].InField = 1;\n')
        f.write(f'Field[2].SizeMin = {refine_size};\n')
        f.write(f'Field[2].SizeMax = {boundary_size};\n')
        f.write('Field[2].DistMin = 0.5;\n')
        f.write('Field[2].DistMax = 5.0;\n\n')
        f.write('Background Field = 2;\n')

    print(f"GEO-скрипт успешно сохранен в файл: {output_filename}")

if __name__ == "__main__":
    CSV_FILE = 'Default Dataset.csv'
    GEO_FILE = 'Mosco-reg.geo'

    BOUNDARY_SIZE = 5.0
    REFINE_SIZE = 0.5

    SHATURA_X = 242.46
    SHATURA_Y = 94.655

    generate_geo_script(CSV_FILE, GEO_FILE, SHATURA_X, SHATURA_Y, REFINE_SIZE, BOUNDARY_SIZE)