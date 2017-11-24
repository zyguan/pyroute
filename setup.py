from setuptools import setup, Extension

setup(name='route-util',
      version='0.0.1',
      packages=['route'],
      author='zyguan',
      license='MIT',
      ext_modules=[Extension('_route', sources=['route.c'],
                             extra_compile_args=['-std=c99'])],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Operating System :: POSIX',
        'Programming Language :: C',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Software Development :: Libraries :: Python Modules',
      ])
