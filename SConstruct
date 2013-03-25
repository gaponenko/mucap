#!/usr/bin/env python
#
# Build a Mu2e base release or test release.
#
# $Id$
# $Author$
# $Date$
#
# Original author Rob Kutschke.
#
import os, re, string
import sys

# Check that the release-specific setup has been run.
if not os.environ.has_key('MU2E_BASE_RELEASE'):
    sys.exit('You must setup a Mu2e base release before running scons.\nExiting.')

# Tell scons about a new command line option that controls the selection of compiler and linker flags.
AddOption('--mu2elevel',
          dest='mu2elevel',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          default='prof',
          help='Select debug build')

# Extract information from the shell environment.
art_inc       = os.environ['ART_INC']
art_lib       = os.environ['ART_LIB']
base          = os.environ['MU2E_BASE_RELEASE']
boost_lib     = os.environ['BOOST_LIB']
boost_inc     = os.environ['BOOST_INC']
clhep_inc     = os.environ['CLHEP_INC']
clhep_lib     = os.environ['CLHEP_LIB_DIR']
cppunit_dir   = os.environ['CPPUNIT_DIR']
gccxml_dir    = os.environ['GCCXML_DIR']
heppdt_lib    = os.environ['HEPPDT_LIB']
heppdt_inc    = os.environ['HEPPDT_INC']
libsigcpp_inc = os.environ['LIBSIGCPP_INC']
libsigcpp_lib = os.environ['LIBSIGCPP_LIB']
root_inc      = os.environ['ROOT_INC']
root_sys      = os.environ['ROOTSYS']
fhicl_inc     = os.environ['FHICLCPP_INC']
fhicl_lib     = os.environ['FHICLCPP_LIB']
cpp0x_inc     = os.environ['CPP0X_INC']
cpp0x_lib     = os.environ['CPP0X_LIB']
mesfac_inc     = os.environ['MESSAGEFACILITY_INC']
mesfac_lib     = os.environ['MESSAGEFACILITY_LIB']
cetlib_inc     = os.environ['CETLIB_INC']
cetlib_lib     = os.environ['CETLIB_LIB']
xercesc_inc    = os.environ['XERCES_C_INC']
xercesc_root   = os.environ['XERCESCROOT']

# If we are working in a test release, extract more information from the environment.
if os.environ.has_key('MU2E_TEST_RELEASE'):
    testrelease          = os.environ['MU2E_TEST_RELEASE']
    cpppath_frag         = [ testrelease, testrelease + '/BaBar/include' ]
    libpath_frag         = [ testrelease+'/lib/' ]
else:
    cpppath_frag         = [ ]
    libpath_frag         = [ ]

# The link libraries needed when building the BaBar code.
babarlibs = [ 'mu2e_BaBar_KalmanTrack',     'mu2e_BaBar_DetectorModel',      'mu2e_BaBar_TrkBase',    'mu2e_BaBar_BField',
              'mu2e_BaBar_TrajGeom',        'mu2e_BaBar_BbrGeom',            'mu2e_BaBar_difAlgebra', 'mu2e_BaBar_ProbTools',
              'mu2e_BaBar_BaBar',           'mu2e_BaBar_CLHEP_src_Geometry', 'mu2e_BaBar_MatEnv',
              'mu2e_BaBar_Dch_DchGeomBase', 'mu2e_BaBar_Dch_DchGeom' ]

# Define scons-local environment - it will be exported later.
osenv = {}
for var in [ 'LD_LIBRARY_PATH',  'GCC_FQ_DIR',  'PATH', 'PYTHONPATH',  'ROOTSYS' ]:
    if var in os.environ.keys():
        osenv[var] = os.environ[var]
        pass
    pass

env = Environment( CPPPATH=[ cpppath_frag,
                             base,
                             base+'/BaBar/include',
                             art_inc,
                             mesfac_inc,
                             fhicl_inc,
                             cetlib_inc,
                             cpp0x_inc,
                             boost_inc,
                             clhep_inc,
                             cppunit_dir+'/include',
                             heppdt_inc,
                             libsigcpp_inc+'/sigc++-2.0',
                             libsigcpp_lib+'/sigc++-2.0/include',
                             root_inc,
                             xercesc_inc,
                           ],
                   LIBPATH=[ libpath_frag,
                             base+'/lib',
                             art_lib,
                             mesfac_lib,
                             fhicl_lib,
                             cetlib_lib,
                             cpp0x_lib,
                             boost_lib,
                             clhep_lib,
                             cppunit_dir+'/lib',
                             heppdt_lib,
                             libsigcpp_lib,
                             root_sys+'/lib',
                             '/lib', '/usr/X11R6/lib',
                             xercesc_root+'/lib',
                           ],
                   ENV=osenv,
                   FORTRAN = 'gfortran',
                   BABARLIBS = [ babarlibs ]
                 )

# Define the rule for building dictionaries.
genreflex_flags = '--deep --fail_on_warnings --iocomments --capabilities=classes_ids.cc '\
                + '-D_REENTRANT -DGNU_SOURCE -DGNU_GCC -D__STRICT_ANSI__ '\
                + '-DPROJECT_NAME="mu2e" -DPROJECT_VERSION="development"'
aa="if   t1=`expr ${TARGET} : '\(.*\)_dict.cpp'`;then t2=$${t1}_map.cpp; t1=$${t1}_dict.cpp;"\
  +"elif t1=`expr ${TARGET} : '\(.*\)_map.cpp'`;then t2=$${t1}_map.cpp; t1=$${t1}_dict.cpp; fi;"\
  +"if genreflex $SOURCE -s ${SOURCE.srcdir}/classes_def.xml $_CPPINCFLAGS"\
  +" -o $$t1 "\
  +genreflex_flags\
  +"; then mv ${TARGET.dir}/classes_ids.cc $$t2; else rm -f $$t1; false; fi"

genreflex = Builder(action=aa)
env.Append(BUILDERS = {'DictionarySource' : genreflex})

# Get the flag that controls compiler options. Check that it is legal.
# There is probably a way to tell AddOption to do this test internally.
level = GetOption('mu2elevel')
known_levels = ['prof', 'debug' ]
if not level in known_levels:
    print 'Unrecognized value for --mu2elevel ' + level
    print '   The value must be one of the known levels: '  + str(known_levels)
    raise Exception('foo')

# Set compile and link flags.
SetOption('warn', 'no-fortran-cxx-mix')
env.MergeFlags('-std=c++11')
env.MergeFlags('-rdynamic')
env.MergeFlags('-Wall')
env.MergeFlags('-g')
if level == 'prof':
    env.MergeFlags('-O3')
    env.MergeFlags('-fno-omit-frame-pointer')
    env.MergeFlags('-DNDEBUG')

if level == 'debug':
    env.MergeFlags('-O0')

# Extract gcc version.  Some libraries have this version embedded in their names.
ff = os.popen('g++ --version'); ll = ff.readline(); ff.close()
gcc_version = ll[10:13]
env.gcc_ver=gcc_version.replace('.','')

# This comes from: root-config --cflags --glibs
# Then guess at the correct location of Spectrum and MLP.
rootlibs = [ 'Core', 'Cint', 'RIO', 'Net', 'Hist', 'Spectrum', 'MLP', 'Graf', 'Graf3d', 'Gpad', 'Tree',
             'Rint', 'Postscript', 'Matrix', 'Physics', 'MathCore', 'Thread', 'Gui', 'm', 'dl' ]
env.Append( ROOTLIBS = rootlibs );

# Make the modified environment visible to all of the SConscript files
Export('env')

# Walk the directory tree to locate all SConscript files.
ss=[]
for root,dirs,files in os.walk('.'):
    for file in files:
        if file == 'SConscript': ss.append('%s/%s'%(root[2:],file))
        pass
    pass

# Define a helper class to construct names of .so libaries. Make an instance of it available to the SConscript files.
class mu2e_helper:
    """mu2e_helper: class to produce library names"""
#   This appears to behave like c++ static member and is initialized at class defintion time.
    sourceroot =  os.path.abspath('.')
#
#   Accesor
#
    def base(self):
        return self.sourceroot
#
#   Build the name of the shared library into which non-plugin compiled code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def libname(self):
        relpath = os.path.relpath('.',self.sourceroot)
        tokens = string.split(relpath,'/')
        if len(tokens) > 1:
            if tokens[len(tokens)-1] == 'src':
                tokens.pop()
                pass
            pass
        return 'mu2e_' + string.join(tokens,'_')
    def prefixed_libname(self):
        return '#/lib/' + self.libname()
#
#   Build the name of the shared library into which plugin code will be inserted.
#   Two versions: with and without the '#/lib' path prefix.
#
    def plugin_libname(self,sourcename):
        return self.libname() + '_' + sourcename[:sourcename.find('.cc')]
    def prefixed_plugin_libname(self,sourcename):
        return '#/lib/' + self.plugin_libname(sourcename)
#
#   Build a list of plugins to be biult.
#
    def plugin_cc(self):
        return Glob('*_module.cc', strings=True) + Glob('*_service.cc', strings=True) + Glob('*_source.cc', strings=True)
#
#   Build a list of .cc files that are not plugings; these go into the library named after the directory.
#
    def non_plugin_cc(self):
        tmp = non_plugin_cc = Glob('*.cc', strings=True)
        for cc in self.plugin_cc(): tmp.remove(cc)
        return tmp
#
#   Names need to build the _dict and _map libraries.
#
    def dict_tmp_name(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return '#/tmp/src/' + relpath + '/' + self.libname() + '_dict.cpp'

    def map_tmp_name(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return '#/tmp/src/' + relpath + '/' + self.libname() + '_map.cpp'

    def dict_libname(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return self.libname() + '_dict'

    def map_libname(self):
        relpath = os.path.relpath('.',self.sourceroot)
        return self.libname() + '_map'

    def prefixed_dict_libname(self):
        return '#/lib/' + self.dict_libname()

    def prefixed_map_libname(self):
        return '#/lib/' + self.map_libname()
#
#   Make the main library.
#
    def make_mainlib( self, userlibs, cppf=[], pf=[], addfortran=False ):
        non_plugin_cc = self.non_plugin_cc()
        if addfortran:
            fortran = Glob('*.f', strings=True)
            non_plugin_cc = [ non_plugin_cc, fortran]
            pass
        libs = []
        if non_plugin_cc:
            env.SharedLibrary( self.prefixed_libname(),
                               non_plugin_cc,
                               LIBS=[ userlibs ],
                               CPPFLAGS=cppf,
                               parse_flags=pf
                              )
            libs = [ self.libname() ]
            pass
        return libs
#
#   Make one plugin library ( but does not work for _dict and _map plugins )
#
    def make_plugin( self, cc, userlibs, cppf = [], pf = []):
        env.SharedLibrary( self.prefixed_plugin_libname(cc),
                           cc,
                           LIBS=[ userlibs, ],
                           CPPFLAGS=cppf,
                           parse_flags=pf
                           )
#
#   Make all plugin libraries, excluding _dict and _map; this works if all libraries need the same link list.
#
    def make_plugins( self, userlibs, exclude_cc = [], cppf = [], pf = [] ):
        plugin_cc = self.plugin_cc()
        for cc in exclude_cc: plugin_cc.remove(cc)
        for cc in plugin_cc:
            env.SharedLibrary( self.prefixed_plugin_libname(cc),
                               cc,
                               LIBS=[ userlibs ],
                               CPPFLAGS=cppf,
                               parse_flags=pf
                               )

#
#   Make the dictionary and map plugins.
#
    def make_dict_and_map( self, userlibs ):
        if os.path.exists('classes.h'):
            if os.path.exists('classes_def.xml'):
                env.DictionarySource([ self.dict_tmp_name(),
                                       self.map_tmp_name() ],
                                     [ 'classes.h', 'classes_def.xml'] )
                env.SharedLibrary( self.prefixed_dict_libname(),
                                   self.dict_tmp_name(),
                                   LIBS=[ userlibs ]
                                   )
                env.SharedLibrary( self.prefixed_map_libname(),
                                   self.map_tmp_name()
                                   )

# Export the class so that it can be used in the SConscript files
# For reasons I don't understand, this must come before the env.SConscript(ss) line.
Export('mu2e_helper')

# Tell scons to operate on all of the SConscript files found by walking the directory tree.
env.SConscript(ss)

# This tells emacs to view this file in python mode.
# Local Variables:
# mode:python
# End:
