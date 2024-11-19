def read_between_lines(file, start, end):
    file.seek(0)
    lines = []
    for i, line in enumerate(file):
        if start <= i < end:
            lines.append(line)
    return "".join(lines)