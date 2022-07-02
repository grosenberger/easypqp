import os

# Global variable for location of packaged unimod.xml file
location = os.path.dirname(os.path.realpath(__file__))
pkg_unimod_db = os.path.join(location, 'data', 'unimod.xml')