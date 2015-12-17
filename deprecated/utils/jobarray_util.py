#!/usr/bin/env python

import os
import sys
import re
import shutil
from optparse import OptionParser

from jobarray_config import config

class jobarray:
    cmd_prefix = 'command_'

    def get_testset_dirs(self, testset_name):
        'Return list of testcase directories for a testset'
        search_dir = config.tc_dir + '/' + testset_name
        findCmd = 'find ' + search_dir + ' -name "' + config.run_dir_search_glob + '" -maxdepth 4 | sort'
        
        ff = os.popen(findCmd, 'r')
        runDirs = [x.strip() for x in ff.readlines()]
        ff.close()
        
        return runDirs

    def get_testcase_list(self, testcase_list_file):
        'Returns testcase list from a file'
        
        testcaseListObj = open(testcase_list_file, 'r')
        runDirs = [x.strip() for x in testcaseListObj.readlines()]
        testcaseListObj.close()
        
        return runDirs

    def get_commands(self, as_names=False):
        'Return list of command functions'
        command_names = []
        for dirName in dir(self):
            if dirName.find(self.cmd_prefix, 0, len(self.cmd_prefix)) != -1:
                if as_names:
                    command_names.append(dirName.replace(self.cmd_prefix, ''))
                else:
                    command_names.append(dirName)            

        return command_names
   
    def print_testset_names(self):
        print 'Valid test set names:'
        print '---------------------'
        for file in os.listdir(config.tc_dir):
            print file

    def process_commands(self):
        # Set up command line arguments
        parser = OptionParser()
           
        parser.add_option( "-l", "--list", dest="list",
                           action="store_true",
                           help="List test case directories for a test set",
                           default=False
                           )

        parser.add_option( "-s", "--script", dest="script",
                           action="store_true",
                           help="Output an LSF script for a particular test set directory or test case list",
                           default=False
                           )
        
        parser.add_option( "-e", "--exec", dest="execute",
                           action="store_true",
                           help="Execute a particular test case for a test case list, binary and index",
                           default=False
                           )

        parser.add_option( "-t", "--testset_name", dest="testset_name",
                           metavar="NAME",
                           help="Name of testset to use for listing or script creation"
                           )

        parser.add_option( "-c", "--testcase_list", dest="testcase_list",
                           metavar="FILE",
                           help="File with list of test case directories. Overrides any list derived soley from testset name"
                           )

        parser.add_option( "-x", "--index", dest="index",
                           metavar="NUM", type="int",
                           help="Index number to either list or execute"
                           )

        parser.add_option( "-b", "--binary", dest="binary_name",
                           metavar="FILE",
                           default=config.binary_name,
                           help="Binary executable to use for testing [default: %default]"
                           )

        (options, args) = parser.parse_args()

        if (options.list and options.script) or \
           (options.script and options.execute) or \
           (options.execute and options.list):
            parser.error('--list, --script and --exec can not be specified together')

        testcase_dirs = []
        if options.testcase_list != None:
            testcase_dirs = self.get_testcase_list(options.testcase_list)
        elif options.testset_name != None:
            testcase_dirs = self.get_testset_dirs(options.testset_name)
        else:
            self.print_testset_names()
            sys.exit(1)

        if options.index != None:
            options.index = options.index - 1
            
            if options.index < 0:
                parser.error('index values start at 1 not 0')
                
            elif options.index >= len(testcase_dirs):
                parser.error('index exceeds size of testcase list')

        if options.list:           
            self.list( testcase_dirs, options.index )
            
        elif options.script:
            
            if options.testset_name == None:
                parser.error('test set name needed for script option')
                
            self.script( options.testset_name, options.testcase_list, testcase_dirs, options.binary_name )
            
        elif options.execute:
            if options.index == None:
                parser.error('index needed to execute script')

            if options.testcase_list == None:
                parser.error('test case list needed for script option')
                
            self.execute(options.testset_name, options.binary_name, options.index, testcase_dirs )
        else:
            parser.error('Must specify one of --list, --script or --exec')

    def list(self, testcase_dirs, index):
        'List all run directories for a testset'

        if index != None:
            print testcase_dirs[index]
        else:
            for test_run_dir in testcase_dirs:
                print test_run_dir

    def script(self, testset_name, testcase_list, testcase_dirs, binary_name):
        'Output a script for processing a testset'

        num_testcases = len(testcase_dirs)

        if sys.argv[0].find('./') == 0:
            script_name = os.getcwd() + sys.argv[0].replace('./', '/')
        elif sys.argv[0].find('/') == 0:
            script_name = sys.argv[0]
        else:
            script_name = os.getcwd() + '/' + sys.argv[0]
                
        print '#!/bin/sh'
        print '#BSUB -n 1'
        print '#BSUB -W ' + str(config.max_run_time)
        print '#BSUB -o "' + config.output_dir + '/stdout.' + testset_name + '.%I"'
        print '#BSUB -e "' + config.output_dir + '/stderr.' + testset_name + '.%I"'
        print '#BSUB -q ' + config.queue_name
        print '#BSUB -J "%s[1-%d]%%%d"' % (testset_name, num_testcases, config.num_simul_jobs)
        print 'date'

        print script_name, '--exec',
        print '--testset_name', testset_name,
        if testcase_list != None:
            print '--testcase_list', testcase_list,
        print '--binary', binary_name, '--index', '$LSB_JOBINDEX'
        
        print 'date'

    def execute(self, testset_name, binary_name, index, testcase_dirs ):
        'Run the L2 binary for a certain indexed run directory in a testset'

        tcNameResult = re.findall(config.testcase_dir_name_re, testcase_dirs[index])
        tcPathResult = re.findall('.*/' + config.testcase_dir_name_re, testcase_dirs[index])
        if len(tcPathResult) <= 0 or len(tcNameResult) <= 0:
            print >>sys.stderr, 'Could not extract testcase name or path'
            sys.exit(0)

        if not os.path.exists(config.local_dir):
            self.output_message('Creating local directory: ' + config.local_dir)
            os.makedirs(config.local_dir)

        msgStr = 'Clearing local directory: ' + config.local_dir
        self.output_message(msgStr)
        os.system(config.rm_bin + " " + config.local_dir + "/*")

        tcName = tcNameResult[0]
        tcOrigDir = tcPathResult[0]
        tcDestDir = config.local_dir + '/' + tcName
        tcRunDir = testcase_dirs[index].replace(tcOrigDir, tcDestDir)

        msgStr = 'Testcase directory: ' + tcOrigDir
        self.output_message(msgStr)
       
        msgStr = 'Copying to local disk: ' + tcOrigDir + ' => ' + tcDestDir
        self.output_message(msgStr)
        cpCmd = config.cp_bin + " -a " + tcOrigDir + " " + tcDestDir
        os.system(cpCmd)

        msgStr = 'Changing directory to: ' + tcRunDir
        self.output_message(msgStr)
        os.chdir( tcRunDir )

        msgStr = 'Executing: ' + binary_name
        self.output_message(msgStr)
        sys.stdout.flush()
        sys.stderr.flush()
        os.system(binary_name)

        if os.path.exists(tcOrigDir):
            msgStr = 'Removing old networked testcase copy: ' + tcOrigDir
            self.output_message(msgStr)
            rmCmd = config.rm_bin + " -rf " + tcOrigDir
            os.system(rmCmd)

        msgStr = 'Copying from local disk: ' + tcDestDir + ' => ' + tcOrigDir
        self.output_message(msgStr)
        cpCmd = config.cp_bin + " -a " + tcDestDir + " " + tcOrigDir
        os.system(cpCmd)

    def output_message(self, msgStr):
        print >>sys.stdout, msgStr
        print >>sys.stdout, ''
        print >>sys.stderr, msgStr
        print >>sys.stderr, ''


# Main
jobarray().process_commands()

