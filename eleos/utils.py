import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import re



# NumPy extension functions

def find_nearest(array, value):
    """Get the closest item to a given value in a Numpy array
    
    Args:
        array (np.ndarray): The array to search in
        value (float): The value to find
        
    Returns:
        int: The index of the closest item in the array
        float: The closest item in the array"""
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def nanaverage(data, weights, axis=None):
    ma = np.ma.MaskedArray(data, mask=np.isnan(data))
    return np.ma.average(ma, weights=weights)


# String / file operation functions

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
    lines = str(string).split("\n")
    ind = " "*spaces_per_level*level
    return "\n".join([ind+l for l in lines])


def write_nums(file, *nums, sep="    ", dp=6):
    for i, n in enumerate(nums):
        file.write(f"{n:.{dp}e}")
        if i != len(nums)-1:
            file.write(sep)
        else:
            file.write("\n")


def generate_ascii_table(title, headers, row_headers, rows):
    headers = [title] + headers
    num_columns = len(headers)

    # Calculate column widths (based on *stripped* text)
    data_col_widths = [
        max(len(remove_ansi_colors(row[i])) for row in rows)
        for i in range(len(headers) - 1)
    ]
    # Insert row-header column
    data_col_widths.insert(
        0, max(len(remove_ansi_colors(h)) for h in row_headers)
    )
    # Ensure column is at least as wide as the header text
    data_col_widths = [
        max(data_col_widths[i], len(remove_ansi_colors(headers[i])))
        for i in range(num_columns)
    ]

    # Box drawing
    top_left_corner = "┌"
    top_right_corner = "┐"
    bottom_left_corner = "└"
    bottom_right_corner = "┘"
    horizontal = "─"
    vertical = "│"
    cross_top = "┬"
    cross_middle = "┼"
    cross_bottom = "┴"

    def pad_cell(cell_text, col_index):
        """
        1) Measure visible length by stripping ANSI codes.
        2) Add trailing spaces to match `data_col_widths[col_index]`.
        3) If the cell is wrapped in a single color code at the start
           and a reset code at the end, include the spaces *inside* that.
        """
        stripped = remove_ansi_colors(cell_text)
        pad_length = data_col_widths[col_index] - len(stripped)

        # If there's no padding needed or no color code, just do normal pad:
        if pad_length <= 0:
            return cell_text

        # Regex to see if the cell starts with some ANSI color code
        # and ends with a reset code.  This is a *naive* approach:
        #   - Looks for something like "\x1b[...m" at start
        #   - Looks for "\x1b[0m" at end
        # If both exist, we insert spaces before the final reset.
        # Otherwise, just append spaces normally.
        pattern_start = re.compile(r'^(?:\x1B\[.*?m)(.*)', re.DOTALL)
        pattern_end = re.compile(r'(.*)(?:\x1B\[0m)$', re.DOTALL)

        match_start = pattern_start.match(cell_text)
        match_end = pattern_end.match(cell_text)

        if match_start and match_end:
            # The cell starts with some color code and ends with a reset code
            # We'll keep them, and insert spaces *inside* those codes.
            # e.g. "\033[31mTEXT\033[0m" => "\033[31mTEXT    \033[0m"
            start_color = cell_text[: match_start.start(1)]  # up to first text
            end_reset = cell_text[match_end.end(1) :]        # from last text to end
            middle_text = match_start.group(1)
            # Now add spaces
            middle_text += " " * pad_length
            # Reconstruct
            return f"{start_color}{middle_text}{end_reset}"
        else:
            # Not a simple wrap; just add trailing spaces as normal
            return cell_text + (" " * pad_length)

    # Helpers for drawing
    def make_row(cells, left, middle, right):
        padded_cells = [pad_cell(str(c), i) for i, c in enumerate(cells)]
        return f"{left} " + f" {middle} ".join(padded_cells) + f" {right}"

    def make_separator(left, middle, right):
        return (
            f"{left}"
            + f"{middle}".join(horizontal * (data_col_widths[i] + 2)
                               for i in range(num_columns))
            + f"{right}"
        )

    # Build table
    table = []
    table.append(make_separator(top_left_corner, cross_top, top_right_corner))
    table.append(make_row(headers, vertical, vertical, vertical))
    table.append(make_separator("├", cross_middle, "┤"))

    for row_header, row in zip(row_headers, rows):
        table.append(make_row([row_header] + row, vertical, vertical, vertical))

    table.append(make_separator(bottom_left_corner, cross_bottom, bottom_right_corner))
    return "\n".join(table)


def get_floats_from_string(string):
    pattern = r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?"
    floats = re.findall(pattern, string)
    return [float(x) for x in floats]


def get_ints_from_string(string):
    pattern = R"[-+]?(?:\d*\.*\d+)"
    floats = re.findall(pattern, string)
    return [int(x) for x in floats]


def format_decimal_hours(decimal_hours):
    hours = int(decimal_hours)
    minutes = (decimal_hours*60) % 60
    seconds = (decimal_hours*3600) % 60
    return "%dh %02dm %02ds" % (hours, minutes, seconds)


def remove_ansi_colors(string):
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    out = ansi_escape.sub('', str(string))
    return out

