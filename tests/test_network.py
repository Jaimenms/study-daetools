import os

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
modelsdir = os.path.join(parentdir,'models')
sys.path.insert(0,modelsdir)


