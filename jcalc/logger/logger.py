import logging


def log_config(log_file):
    """ Description:
          Function to start log, writing inside LOGS

        Usage:
          log_config("/home/ubuntu/LOGS/jcalc.log")

        Parameters:
          log_file:
            string, LOGS directory from jcalc
    """
    # Log configurations
    logging.basicConfig(
        filename=log_file, filemode="w",
        format="%(asctime)s - %(message)s",
        datefmt='%m/%d/%Y %H:%M:%S',
        level=logging.INFO
    )
