import time
import logging
import functools


def log_function(function):
    func_name = function.__name__
    logger = logging.getLogger()

    @functools.wraps(function)
    def log_function_wrapper(*args, **kwargs):
        logger.info(f">>>>> ENTERING the function '{func_name}(...)'")
        result = function(*args, **kwargs)
        logger.info(f"<<<<< EXITING the function '{func_name}(...)'")
        return result
    
    return log_function_wrapper


class ExecTimeCM(object):
    """Context manager to be used with the 'with' clause 
    to get the execution time of a code block.
    
    Sample Usage:
    .. code:: python
        with ExecTimeCM("Testing") as st:
            print("sample code block.")
    .. code:: console
        
        EXECUTION TIME for Testing: 0.0001 seconds
    """
    
    def __init__(self, txt="", verbose=True):
        """Called when the context manager scope is initiated.
        """
        self.txt = txt
        self.verbose = verbose
        self.description = f"EXECUTION TIME for {self.txt}"
    
    def __enter__(self):
        """Called when the scope is entered.
        """
        self.start = time.time()
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        """Called when the scope is exited.
        """
        self.exec_time = time.time() - self.start
        if self.verbose:
            print(f"{self.description}: {self.exec_time} seconds")