
# system modules

# custom modules
import libpsf.utils.m_timers as m_timers
import libpsf.utils.m_logging as m_logging
import libpsf.input.m_config as m_config

# --------------------------------------------------------------------------------------------------
# set everything up
# --------------------------------------------------------------------------------------------------

# logger to handle printing (and crashing ...)
logger = m_logging.c_logger()
logger.print_herald()

# get options from cmd line
config = m_config.c_config(logger)

# --------------------------------------------------------------------------------------------------
# loop over the files
# --------------------------------------------------------------------------------------------------

# loop over the input files
for ii in range(config.num_files):
    
    input_file = config.input_files[ii]
    msg = f'now doing file \'{input_file}\''
    logger.console(msg)

    # get additional config info from file
    config.get_args(input_file)

# --------------------------------------------------------------------------------------------------
# finalize everything
# --------------------------------------------------------------------------------------------------

# clean everything up
logger.clean_up()



