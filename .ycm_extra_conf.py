import os

def FlagsForFile(filename, **kwargs):

    flags = ['-std=c++14', '-I/usr/local/include']

    proj_root = os.path.dirname(os.path.abspath(__file__))
    proj_include = ''.join(['-I', proj_root, '/include'])
    flags.append(proj_include)

    if filename.endswith('_test.cpp'):
        gtest_include = ''.join([
            '-I', proj_root, '/test/googletest/googletest/include'
            ])
        flags.append(gtest_include)

    return {'flags': flags}
