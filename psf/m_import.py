
# system modules
import importlib
import importlib.machinery as machinery

# custom modules
from psf.m_error import crash


# --------------------------------------------------------------------------------------------------

def import_module(module_name):

    """
    import a module that is located in the source 'tree' 
    """

    try:
        module = importlib.import_module(module_name)
    except Exception as _ex:
        crash(f'the file \'{module_name}\' is broken\n',_ex)

    return module

# --------------------------------------------------------------------------------------------------

def import_module_from_path(module_path):

    """
    import a module with arbitrary path. hopefully this method will be robust...
    """

    try:
        module = importlib.machinery.SourceFileLoader('module',module_path)
        module = module.load_module()
    except Exception as _ex:
        crash(f'the file \'{module_path}\' is broken\n',_ex)

    return module

# --------------------------------------------------------------------------------------------------




