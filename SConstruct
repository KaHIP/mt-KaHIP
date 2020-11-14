#/******************************************************************************
# * SConstruct
# *
# * Source of KaHIP -- Karlsruhe High Quality Partitioning.
# *
# ******************************************************************************
# * Copyright (C) 2015 Christian Schulz <christian.schulz@kit.edu>
# *
# * This program is free software: you can redistribute it and/or modify it
# * under the terms of the GNU General Public License as published by the Free
# * Software Foundation, either version 5 of the License, or (at your option)
# * any later version.
# *
# * This program is distributed in the hope that it will be useful, but WITHOUT
# * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# * more details.
# *
# * You should have received a copy of the GNU General Public License along with
# * this program.  If not, see <http://www.gnu.org/licenses/>.
# *****************************************************************************/



# scons build file for the KaHIP.
#
# You can build it in the following variants:
#
#   optimized            no debug symbols, no assertions, optimization.
#   optimized_output     no debug symbols, no assertions, optimization -- more output on console.
#
#   scons variant=${variant} program=${program}
import os
import platform
import sys

# Get the current platform.
SYSTEM = platform.uname()[0]
HOST = platform.uname()[1]

# Get shortcut to $HOME.
HOME = os.environ['HOME']

def GetEnvironment():
  """Get environment variables from command line and environment.

  Exits on errors.

  Returns
    Environment with the configuration from the command line.
  """
  opts = Variables()
  opts.Add('variant', 'the variant to build, optimized or optimized with output', 'optimized')
  opts.Add('program', 'program or interface to compile', 'kaffpa')

  env = Environment(options=opts, ENV=os.environ)
  if not env['variant'] in ['optimized','optimized_output','debug']:
    print 'Illegal value for variant: %s' % env['variant']
    sys.exit(1)
  
  if not env['program'] in ['kaffpa', 'kaffpa_test', 'kaffpa_compare_with_sequential', 'kaffpa_test_stopping_rule', 'kaffpaE', 'partition_to_vertex_separator','improve_vertex_separator','library','graphchecker','label_propagation','evaluator','node_separator']:
    print 'Illegal value for program: %s' % env['program']
    sys.exit(1)

  # Special configuration for 64 bit machines.
  if platform.architecture()[0] == '64bit':
     env.Append(CPPFLAGS=['-DPOINTER64=1'])

  return env

# Get the common environment.
env = GetEnvironment()

env.Append(LIBPATH=['/home/akhremtsev/lib/tbb/build/linux_intel64_gcc_cc5.4.1_libc2.19_kernel3.13.0_release'])
env.Append(CPPPATH=['../extern/argtable-2.10/include'])
env.Append(CPPPATH=['/opt/intel/vtune_amplifier_2018/include'])
env.Append(CPPPATH=['./extern/argtable-2.10/include'])
env.Append(CPPPATH=['./lib'])
env.Append(CPPPATH=['/home/akhremtsev/lib/tbb/include'])
env.Append(CPPPATH=['./lib/tools'])
env.Append(CPPPATH=['./lib/partition'])
env.Append(CPPPATH=['./lib/io'])
env.Append(CPPPATH=['./lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/'])
env.Append(CPPPATH=['./growt/'])
env.Append(CPPPATH=['./ips4o/'])
env.Append(CPPPATH=['./stats/include/'])
env.Append(LIBPATH=['../extern/argtable-2.10/lib'])
env.Append(LIBPATH=['./extern/argtable-2.10/lib'])
env.Append(CPPPATH=['../lib'])
env.Append(CPPPATH=['../lib/tools'])
env.Append(CPPPATH=['../lib/partition'])
env.Append(CPPPATH=['../lib/io'])
env.Append(CPPPATH=['../lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/'])
env.Append(LIBPATH=['../../extern/argtable-2.10/lib'])
env.Append(LIBPATH=['/opt/intel/vtune_amplifier_2018/lib64/'])
env.Append(CPPPATH=['/usr/include/openmpi/'])

conf = Configure(env)

if SYSTEM == 'Darwin':
        env.Append(CPPPATH=['/opt/local/include/','../include'])
        env.Append(LIBPATH=['/opt/local/lib/'])
        env.Append(LIBPATH=['/opt/local/lib/openmpi/'])
        # homebrew related paths
        env.Append(LIBPATH=['/usr/local/lib/'])
        env.Append(LIBPATH=['/usr/local/lib/openmpi/'])

#by D. Luxen
#if not conf.CheckLibWithHeader('argtable2', 'argtable2.h', 'CXX'):
        #print "argtable library not found. Exiting"
        #Exit(-1)
#if not conf.CheckCXXHeader('mpi.h'):
        #print "openmpi header not found. Exiting"
        #Exit(-1)
#
#
env.Append(CXXFLAGS = '-fopenmp -g')
#env.Append(CXXFLAGS = '-g')
# Apply variant specific settings.
if env['variant'] == 'optimized':
  env.Append(CXXFLAGS = '-DNDEBUG -Wall -funroll-loops  -fno-stack-limit -O3 -std=c++1z -g -mcx16')
  env.Append(CCFLAGS  = '-O3  -DNDEBUG -funroll-loops -std=c++1z')
elif env['variant'] == 'optimized_output':
  # A little bit more output on the console
  env.Append(CXXFLAGS = ' -DNDEBUG -funroll-loops -Wall -fno-stack-limit -O3 -std=c++14 -g')
  env.Append(CCFLAGS  = '-O3  -DNDEBUG -DKAFFPAOUTPUT  -std=c++14')
else:
  #env.Append(CXXFLAGS = ' -DNDEBUG -Wall -funroll-loops  -fno-stack-limit -O3 -std=c++14 -g')
  #env.Append(CCFLAGS  = '-O3  -DNDEBUG -funroll-loops -std=c++14 ')
  env.Append(CXXFLAGS = '-O0 -std=c++14 -g')
  if SYSTEM != 'Darwin':
        env.Append(CXXFLAGS = '-march=native')
        env.Append(CCFLAGS  = '-march=native')

# Execute the SConscript.
SConscript('SConscript', exports=['env'],variant_dir=env['variant'], duplicate=False)

