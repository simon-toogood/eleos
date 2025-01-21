def read_between_lines(file, start, end):
    """Take in an open file object and get the data between line numbers 'start' and 'end'.
    It is inclusive of start and exclusive of end (ie. start <= line < end)
    
    Args:
        file: File-like object to read from
        start: Line number to start reading from
        end: Line number to stop reading at
        
    Returns:
        str: File data
    """
    file.seek(0)
    lines = []
    for i, line in enumerate(file):
        if start <= i < end:
            lines.append(line)
    return "".join(lines)


def indent(string, level=1, spaces_per_level=4):
    """Take in a string and return the same string indented at each line break.
    
    Args:
        string: The string to indent
        level: The level at which to indent
        spaces_per_level: How many spaces to add per level of indentation
        
    Returns:
        str: The indented string
    """
    lines = string.split("\n")
    ind = " "*spaces_per_level*level
    return "\n".join([ind+l for l in lines])


def write_nums(file, *nums, sep="    ", dp=6):
    for i, n in enumerate(nums):
        file.write(f"{n:.{dp}f}")
        if i != len(nums)-1:
            file.write(sep)
        else:
            file.write("\n")


def generate_ascii_table(title, headers, row_headers, rows):
    """Create a text-based table with a header column using Unicode box drawing characters.

    Args:
        headers (list): A list of column headers.
        row_headers (list): A list of row header labels (one for each row).
        rows (list of list): A list of rows, where each row is a list of column values.
        top_left (str): The string to place in the top-left cell (default is empty).

    Returns:
        str: The formatted table as a string.
    """
    # Add the top-left string to the headers
    headers = [title] + headers

    # Determine the column widths, including the row header column
    num_columns = len(headers)
    col_widths = [max(len(str(row[i])) for row in rows) for i in range(len(headers) - 1)]
    col_widths.insert(0, max(len(str(h)) for h in row_headers))  # Row header column
    col_widths = [max(col_widths[i], len(headers[i])) for i in range(num_columns)]

    # Box drawing characters
    top_left_corner = "┌"
    top_right_corner = "┐"
    bottom_left_corner = "└"
    bottom_right_corner = "┘"
    horizontal = "─"
    vertical = "│"
    cross_top = "┬"
    cross_middle = "┼"
    cross_bottom = "┴"

    # Helper to create rows
    def make_row(cells, left, middle, right):
        return f"{left} " + f" {middle} ".join(f"{str(cell):<{col_widths[i]}}" for i, cell in enumerate(cells)) + f" {right}"

    # Header separator
    def make_separator(left, middle, right):
        return f"{left}" + f"{middle}".join(horizontal * (col_widths[i] + 2) for i in range(num_columns)) + f"{right}"

    # Build the table
    table = [
        make_separator(top_left_corner, cross_top, top_right_corner),  # Top border
        make_row(headers, vertical, vertical, vertical),  # Header row
        make_separator("├", cross_middle, "┤"),  # Separator under header
    ]
    for row_header, row in zip(row_headers, rows):
        table.append(make_row([row_header] + row, vertical, vertical, vertical))
    table.append(make_separator(bottom_left_corner, cross_bottom, bottom_right_corner))  # Bottom border

    return "\n".join(table)