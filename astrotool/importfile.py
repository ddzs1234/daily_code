"""
multiple ways to import file
"""

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/path/to/application/app/folder')
import file

# https://stackoverflow.com/questions/4383571/importing-files-from-different-folder