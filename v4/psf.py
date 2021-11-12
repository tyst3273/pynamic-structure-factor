
# system modules

# custom modules
import libpsf.utils.m_timers as m_timers
import libpsf.utils.m_logging as m_logging
import libpsf.io.m_config as m_config
import libpsf.structure.m_lattice as m_lattice
import libpsf.structure.m_qpoints as m_qpoints
import libpsf.calc.m_trajectory as m_trajectory
import libpsf.calc.m_xlengths as m_xlengths
import libpsf.calc.m_sqw as m_sqw

# --------------------------------------------------------------------------------------------------
# set everything up
# --------------------------------------------------------------------------------------------------

# set up timer
timers = m_timers.c_timers()

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

    # init lattice
    lattice = m_lattice.c_lattice(config,logger)

    # init Q-points
    qpts = m_qpoints.c_qpoints(config,logger,lattice)

    # init trajectory file
    traj = m_trajectory.c_trajectory(config,logger,lattice)

    # init scattering lenghts
    xlens = m_xlengths.c_xlengths(config,logger)

    # init sqw calculator
    calc = m_sqw.c_sqw(config,logger,lattice,qpts,traj,xlens,timers)

    # run the calculation 
    calc.run_calc()

# --------------------------------------------------------------------------------------------------
# finalize everything
# --------------------------------------------------------------------------------------------------

# print goodbye info
logger.print_goodbye()

# clean everything up
logger.clean_up()




