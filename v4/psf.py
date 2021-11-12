
# system modules

# custom modules
import libpsf.utils.m_timers as m_timers
import libpsf.utils.m_logging as m_logging
import libpsf.input.m_config as m_config

# --------------------------------------------------------------------------------------------------

# logger to handle
logger = m_logging.c_logger()
logger.print_herald()

# get config info from cmd line
config = m_config.c_config(logger)

# --------------------------------------------------------------------------------------------------

# loop over the input files
for ii in range(config.num_files):
    
    input_file = config.input_files[ii]
    msg = f'now doing file \'{input_file}\''
    logger.console(msg)

    config.get_args(input_file)

# --------------------------------------------------------------------------------------------------

# clean everything up
logger.clean_up()



