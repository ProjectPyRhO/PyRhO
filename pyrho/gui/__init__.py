from pkg_resources import resource_filename

# resource_filename(package_or_requirement, resource_name)

# TODO: Boilerplate for widgets and notebook environment
try:
    import ipywidgets as widgets
    from IPython.display import display
    from IPython.display import clear_output
except ImportError:
    widgets = None
