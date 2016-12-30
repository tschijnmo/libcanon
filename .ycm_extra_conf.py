import os

def FlagsForFile(filename, **kwargs):
    proj_root = os.path.dirname(os.path.abspath(__file__))
    proj_include = ''.join(['-I', proj_root, '/include'])

    return {'flags': [
        '-std=c++1z',
        '-I/usr/local/include',
        proj_include
        ]}
