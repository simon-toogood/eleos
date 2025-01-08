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