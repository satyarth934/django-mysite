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