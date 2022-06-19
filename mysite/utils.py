def warning_format(msg, *args, **kwargs):
    """Used to generate warnings within code.
    Prints the warning in the following format:
    <script>:<n>: <type of warning>: <warning msg>

    'n' is the line number of initiated warning command within the script

    Example Usage:
        import warnings
        warnings.formatwarning = warning_format
        warnings.warn("<Warning msg>")
    """

    return f"{args[1]}:{args[2]}: {args[0].__name__}: {str(msg)}\n"