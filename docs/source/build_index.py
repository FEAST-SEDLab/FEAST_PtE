import os
import sys
import importlib
import inspect
sys.path.append('../../')

section_symbols = ['=', '-', '*', '+']
sec_symbol_length = 40

f = open('index.rst', 'w')

lines = ["Welcome to FEAST's documentation!\n"]
lines.append(section_symbols[0] * sec_symbol_length + '\n')
lines.append("\n .. toctree::\n    :maxdepth: 5\n    :caption: Contents:\n \nFEAST modules\n")
lines.append(section_symbols[0] * sec_symbol_length + '\n\n')


def read_module(name, level=1, path='feast.', path_in='../../feast'):
    if '__' in name or ('.' in name and '.py' not in name):
        return
    elif '.py' in name:
        lines.append(name[:-3] + '\n' + section_symbols[level] * sec_symbol_length + '\n')
        lines.append('.. automodule:: ' + path + name[:-3] + '\n\n')
        mod = importlib.import_module(path + name[:-3])
        for cl_name, obj in inspect.getmembers(mod):
            if inspect.isclass(obj):
                # This if condition prevents imported classes from being printed twice
                if obj.__module__ == mod.__name__:
                    lines.append(obj.__name__ + '\n' + section_symbols[level + 1] * sec_symbol_length + '\n')
                    lines.append('.. autoclass:: ' + path + name[:-3] + '.' + obj.__name__ + '\n')
                    lines.append('    :members:\n\n')
            elif inspect.isfunction(obj):
                # This if condition prevents imported classes from being printed twice
                if obj.__module__ == mod.__name__:
                    lines.append(obj.__name__ + '\n' + section_symbols[level + 1] * sec_symbol_length + '\n')
                    lines.append('.. autofunction:: ' + path + name[:-3] + '.' + obj.__name__ + '\n\n')
    else:
        lines.append(name + '\n' + section_symbols[level] * sec_symbol_length + '\n')
        for new_name in os.listdir(os.path.join(path_in, name)):
            read_module(new_name, level=level + 1, path=path + name + '.', path_in=os.path.join(path_in, name))


for name in os.listdir('../../feast'):
    read_module(name)


f.writelines(lines)

f.close()
