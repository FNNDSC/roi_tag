#!/usr/bin/env python

'''

    'roi_controller.py' is a high level "controller" that runs the necessary
    steps to (sub)sample an experimental space and redo the entire cloud
    and tag analysis.

'''


import  os
import  sys
import  string
import  argparse
import  tempfile, shutil
import  json
import  re
import  time

from    _common import systemMisc       as misc
from    _common import crun

import  error
import  message
import  stage
import  fnndsc  as base
import  C_snode


class FNNDSC_ROICONTROLLER(base.FNNDSC):
    '''
    This class is a specialization of the FNNDSC base and geared to dyslexia
    curvature analysis.

    '''

    #
    # Class member variables -- if declared here are shared
    # across all instances of this class
    #
    _dictErr = {
        'subjectSpecFail'   : {
            'action'        : 'examining command line arguments, ',
            'error'         : 'it seems that no subjects were specified.',
            'exitCode'      : 10},
        'noFreeSurferEnv'   : {
            'action'        : 'examining environment, ',
            'error'         : 'it seems that the FreeSurfer environment has not been sourced.',
            'exitCode'      : 11},
        'noStagePostConditions' : {
            'action'        : 'querying a stage for its exitCode, ',
            'error'         : 'it seems that the stage has not been specified.',
            'exitCode'      : 12},
        'subjectDirnotExist': {
            'action'        : 'examining the <subjectDirectories>, ',
            'error'         : 'the directory does not exist.',
            'exitCode'      : 13},
        'Load'              : {
            'action'        : 'attempting to pickle load object, ',
            'error'         : 'a PickleError occured.',
            'exitCode'      : 14},
        'outDirNotCreate': {
            'action'        : 'attempting to create the <outDir>, ',
            'error'         : 'a system error was encountered. Do you have create permission?',
            'exitCode'      : 15},
        'outDirNotFound': {
            'action'        : 'checking on the <outDir>, ',
            'error'         : "directory doesn't seem to exist. Will attempt to create.",
            'exitCode'      : 16},
        'subsampleDirNotFound': {
            'action'        : 'examining a subsample directory, ',
            'error'         : 'a system error was encountered -- directory not found.',
            'exitCode'      : 17},
        'workingDirNotExist': {
            'action'        : 'attempting to access the <workingDir>, ',
            'error'         : 'a system error was encountered. Does the directory exist?',
            'exitCode'      : 18},
    }

    def dtree(self):
        return self._dtree

    def l_pval(self):
        return self._l_pval

    def l_roi(self):
        return self._l_ROI

    def l_subsample(self):
        return self._l_subsample

    def l_hemisphere(self):
        return self._l_hemi

    def l_surface(self):
        return self._l_surface

    def l_statFunc(self):
        return self._l_statFunc

    def l_group(self):
        return self._l_group

    def l_subsample(self):
        return self._l_subsample

    def l_subj(self):
        return self._l_subj

    def l_subsampleDir(self):
        return self._l_subsampleDir

    def l_ctype(self):
        return self._l_ctype

    def l_curvFunc(self):
        return self._l_curvFunc

    def l_annotation(self):
        return self._l_annotation

    def pval(self):
        return self._str_pval

    def topDir(self, *args):
        if len(args):
            self._topDir    = args[0]
        else:
            return self._topDir

    def dtree_build(self):
        '''
        Constructs the internal snode tree for holding data sample points.
        '''

        c = self._dtree
        c.mknode(self._l_annotation)
        for self._str_annotation in self._l_annotation:
            c.cdnode('/')
            c.cdnode(self._str_annotation)
            c.mknode(self._l_group)
            for self._str_group in self._l_group:
                c.cdnode(self._str_group)
                c.mknode(self._l_pval)
                for self._str_pval in self._l_pval:
                    c.cdnode(self._str_pval)
                    c.mknode(self._l_statFunc)
                    for self._str_statFunc in self._l_statFunc:
                        c.cdnode(self._str_statFunc)
                        c.mknode(self._l_surface)
                        for self._str_surface in self._l_surface:
                            c.cdnode(self._str_surface)
                            c.mknode(self._l_hemi)
                            for self._str_hemi in self._l_hemi:
                                c.cdnode(self._str_hemi)
                                c.mknode(self._l_ctype)
                                for self._str_ctype in self._l_ctype:
                                    c.cdnode(self._str_ctype)
                                    c.mknode(self._l_subsample)
                                    c.cdnode('../')
                                c.cdnode('../')
                            c.cdnode('../')
                        c.cdnode('../')
                    c.cdnode('../')
                c.cdnode('../')
            c.cdnode('../')
        c.cdnode('../')

    def dtreeDir(self):
        '''
        Return the 'stree' dir based on internal loop state.
        '''
        str_dtreeDir    = '/%s/%s/%s/%s/%s/%s/%s/%s' % (
                            self._str_annotation,
                            self._str_group,
                            self._str_pval,
                            self._str_statFunc,
                            self._str_surface,
                            self._str_hemi,
                            self._str_ctype,
                            self._str_subsample
            )
        return str_dtreeDir

    def d2path(self):
        '''
        Converts the current tree working 'dir' to a location
        on the hard drive.
        '''
        str_hdPath  = '%s/%s' % (
                            self.outDir(),
                            self._dtree.cwd()
            )
        return str_hdPath

    def dirSpecSubsample(self):
        '''
        Return the directory spec of the current sub-sample
        experiment.
        '''
        str_dirSpecSubsample = '%s/groupCurvAnalysis/%s/ROItag/%s/%s/%s/%s/%s/%s' % (
                                        self._d_subsampleDir[self._str_subsample],
                                        self._str_annotation,
                                        self._str_group,
                                        self._str_pval,
                                        self._str_statFunc,
                                        self._str_surface,
                                        self._str_hemi,
                                        self._str_ctype)
        return str_dirSpecSubsample

    def dirSpec(self):
        '''
        Return the dirSpec based on internal pipeline._str_* variables
        '''
        return '%s/%s/%s/%s/%s/%s/%s' % (  self.outDir(),
                                        self._str_annotation,
                                        self._str_group,
                                        self._str_pval,
                                        self._str_statFunc,
                                        self._str_surface,
                                        self._str_hemi)

    def dirSpecPartial(self):
        '''
        Return the dirSpec based on internal pipeline._str_* variables w/o
        the leading directories.
        '''
        return '%s/%s/%s/%s/%s/%s' % ( self._str_annotation,
                                    self._str_group,
                                    self._str_pval,
                                    self._str_statFunc,
                                    self._str_surface,
                                    self._str_hemi)


    def roi(self):
        return self._str_roi

    def subsample(self):
        return self._str_subsample

    def surface(self):
        return self._str_surface

    def hemi(self):
        return self._str_hemi

    def statFunc(self):
        return self._str_statFunc

    def outDir(self, *args):
        if len(args):
            self._outDir    = args[0]
        else:
            return self._outDir

    def workingDir(self, *args):
        if len(args):
            self._workingDir = args[0]
        else:
            return self._workingDir

    def schedulerStdOutDir(self, *args):
        if len(args):
            self._str_schedulerStdOutDir    = args[0]
        else:
            return self._str_schedulerStdOutDir

    def schedulerStdErrDir(self, *args):
        if len(args):
            self._str_schedulerStdErrDir    = args[0]
        else:
            return self._str_schedulerStdErrDir

    def group(self):
        return self._str_group

    def __init__(self, **kwargs):
        '''
        Basic constructor. Checks on named input args, checks that files
        exist and creates directories.

        '''
        base.FNNDSC.__init__(self, **kwargs)

        self._lw                        = 120
        self._rw                        = 20
        self._l_ROI                     = []
        self._l_subj                    = []
        self._l_subsample               = []
        self._l_subsampleDir            = []
        self._l_annotation              = []

        # Tracking dictionaries
        self._d_subsampleDir            = {}
        self._d_subsampleSet            = {}

        self._l_pval                    = []
        self._l_group                   = []
        self._l_surface                 = []
        self._l_statFunc                = []
        self._l_curvFunc                = []
        self._l_hemi                    = []
        self._l_ctype                   = ['thickness', 'curv']

        self._outDir                    = ''
        self._workingDir                = ''
        self._stageslist                = '12'

        self._str_threshold             = "0.5"

        # Internal tracking vars
        self._str_pval                  = ''
        self._str_group                 = ''
        self._str_roi                   = ''
        self._str_hemi                  = ''
        self._str_surface               = ''
        self._str_statFunc              = ''
        self._str_subsample             = ''
        self._str_subj                  = ''
        self._str_annotation            = ''
        self._str_ctype                 = ''

        self._topDir                    = ''

        # Scheduler std out/err dirs
        self._str_schedulerStdOutDir    = '~/scratch'
        self._str_schedulerStdErrDir    = '~/scratch'

        self._dtree                     = C_snode.C_stree(['/'])

        for key, value in kwargs.iteritems():
            if key == 'subjectList':    self._l_subj            = value
            if key == 'sampleList':     self._l_subsample       = value.split(',')
            if key == 'annot':          self._l_annotation      = value.split(',')
            if key == 'outDir':         self._outDir            = value
            if key == 'workingDir':     self._workingDir        = value
            if key == 'stages':         self._stageslist        = value
            if key == 'pval':           self._l_pval            = value.split(',')
            if key == 'group':          self._l_group           = value.split(',')
            if key == 'surface':        self._l_surface         = value.split(',')
            if key == 'statFunc':       self._l_statFunc        = value.split(',')
            if key == 'curvFunc':       self._l_curvFunc        = value.split(',')
            if key == 'hemi':           self._l_hemi            = value.split(',')
            if key == 'threshold':      self._str_threshold     = value
            if key == 'schedulerStdOutDir':     self._str_schedulerStdOutDir = value
            if key == 'schedulerStdErrDir':     self._str_schedulerStdErrDir = value

        if not os.path.isdir(self._workingDir): errorFatal(self, 'workingDirNotExist')
        self._topDir     = self._workingDir
        self._l_stages   = list(self._stageslist)

    def initialize(self):
        '''
        This method provides some "post-constructor" initialization. It is
        typically called after the constructor and after other class flags
        have been set (or reset).

        '''

        # Set the stages
        self._pipeline.stages_canRun(False)
        for index in self._l_stages:
            stage = self._pipeline.stage_get(int(index))
            stage.canRun(True)
            stage.exitCode(False)

        # Set absolute dir specs for topDir and outDir
        log = self.log()
        log('Setting dir specs...\n')
        os.chdir(self.topDir())
        self.topDir(os.path.abspath(self.topDir()))

    def run(self):
        '''
        The main 'engine' of the class.

        '''
        base.FNNDSC.run(self)


def synopsis(ab_shortOnly = False):
    scriptName = os.path.basename(sys.argv[0])
    shortSynopsis =  '''
    SYNOPSIS

            %s                                \\
                            [-s|--stages <stages>]          \\
                            [-o|--outDir <outputRootDir>]   \\
                            [-w|--workingDir <workingDir>]  \\
                            [-r|--sample <sampleList>]      \\
                            [-t|--threshold <threshold>]    \\
                            [-a|--annot <annotation>]       \\
                            [-v|--verbosity <verboseLevel>] \\
                            [-p|--pval <pvalCutoffList>]    \\
                            [-g|--group <groupList>]        \\
                            [-S|--surface <surfaceList>]    \\
                            [-c|--curvFunc <curvFuncList>]  \\
                            [-f|--statFunc <statFuncList>]  \\
                            [-m|--hemi <hemiList>]          \\
                            [--subShellFlush]               \\
                            [--subShellDontRunCmd]          \\
                            [--schedulerStdOutDir <dir>]    \\
                            [--schedulerStdErrDir <dir>]    \\
                            <subj1> <subj2> ... <subjN>
    ''' % scriptName

    description =  '''
    DESCRIPTION

        `%s' is a "meta" controller for re-running the analysis pipeline
        on targetted (sub)samples of the original data space.

    ARGS

       --stages <stages>
       The stages to execute. This is specified in a string, such as '1234'
       which would imply stages 1, 2, 3, and 4.

       The special keyword 'all' can be used to turn on all stages.

       --outDir <outputRootDir>
       The stem for the output top level directory. If there is a <sampleList>,
       then each element of the <sampleList> is appended to <outputRootDir>
       to create the actual output directory.

       --workingDir <workingDir>
       The working directory for the script.

       --sample <sampleList>
       A comma separated list of sub-sample text suffices.

       --annot <annotation>
       The annotation to analyze. Assumes that this directory contains a subtree,
       'ROItag'.

       --pval <pvalCutoffList>
       The pval cutoffs to consider. In practice, this is always 'le1,le5'

       --group <groupList>
       The group list to process.

       --surface <surfaceList>
       The surface list to process. In practice, 'smoothwm,pial'

       --statFunc <statFuncList>
       The statistical functional data to analyze. Typically
       'ptile-raw,ptile-convex'

       --curvFunc <curvFuncList>
       The curvature functions to ROI tag.

       --hemi <hemiList>
       The hemispheres to process. In practice, this is always 'lh,rh'.

       --subShellFlush
       If passed, flush outputs of internal subshell cmds. Useful mainly for debugging.

       --subShellDontRunCmd
       If passed, do NOT run any subshell cmds. Useful mainly for debugging.

       --schedulerStdOutDir <dir>
       The directory to contain stdout from any scheduled children jobs.

       --schedulerStdErrDir <dir>
       The directory to contain stderr from any scheduled children jobs.

       <subj1> <subj2> ... <subjN>
       The subjects to process.

    EXAMPLES


    ''' % (scriptName)
    if ab_shortOnly:
        return shortSynopsis
    else:
        return shortSynopsis + description

def f_stageShellExitCode(**kwargs):
    '''
    A simple function that returns a conditional based on the
    exitCode of the passed stage object. It assumes global access
    to the <pipeline> object.

    **kwargs:

        obj=<stage>
        The stage to query for exitStatus.

    '''
    stage = None
    for key, val in kwargs.iteritems():
        if key == 'obj':                stage                   = val
    if not stage: error.fatal(pipeline, "noStagePostConditions")
    if not stage.callCount():   return True
    if not stage.exitCode():    return True
    else:                       return False

def l_filterNumeric(l_input):
    '''
    Given an input string list, return a list containing only those original
    elements that have a string numeric component.
    '''
    l_filter = []
    b_add    = False
    for el in l_input:
        try:
            numeric = re.search(r'\d+', el).group()
            b_add = True
        except:
            b_add = False
        if b_add: l_filter.append(el)
    return l_filter

#
# entry point
#
if __name__ == "__main__":


    # always show the help if no arguments were specified
    if len( sys.argv ) == 1:
        print synopsis()
        sys.exit( 1 )

    verbosity   = 0

    parser = argparse.ArgumentParser(description = synopsis(True))

    parser.add_argument('l_subj',
                        metavar='SUBJs', nargs='+',
                        help='Subject directories to process')
    parser.add_argument('--subShellFlush',
                        dest='subShellFlush',
                        action='store_true',
                        default=False,
                        help='flush sub-shell outputs as they run -- mainly debugging')
    parser.add_argument('--subShellDontRunCmd',
                        dest='subShellDontRunCmd',
                        action='store_true',
                        default=False,
                        help='run (or not) all subshell commands -- mainly debugging')
    parser.add_argument('--sample', '-r',
                        dest='sample',
                        action='store',
                        default='',
                        help='comma separated resampled experiment suffix list to process')
    parser.add_argument('--threshold', '-t',
                        dest='threshold',
                        action='store',
                        default='0.50',
                        help='directory selection threshold')
    parser.add_argument('--verbosity', '-v',
                        dest='verbosity',
                        action='store',
                        default=0,
                        help='verbosity level')
    parser.add_argument('--annot', '-a',
                        dest='annotation',
                        action='store',
                        default='',
                        help='comma separated annotation list to process')
    parser.add_argument('--output', '-o',
                        dest='outDir',
                        action='store',
                        default='roicontroller',
                        help='output root directory')
    parser.add_argument('--workingDir', '-w',
                        dest='workingDir',
                        action='store',
                        default='./',
                        help='output working directory')
    parser.add_argument('--stages', '-s',
                        dest='stages',
                        action='store',
                        default='0',
                        help='analysis stages')
    parser.add_argument('--pval', '-p',
                        dest='pval',
                        action='store',
                        default='le1',
                        help='comma separated p-val cutoff threshold')
    parser.add_argument('--group', '-g',
                        dest='group',
                        action='store',
                        default='13',
                        help='comma separated group list to process')
    parser.add_argument('--surface', '-S',
                        dest='surface',
                        action='store',
                        default='smoothwm',
                        help='comma separated surface list to process')
    parser.add_argument('--curvFunc', '-c',
                        dest='curvFunc',
                        action='store',
                        default='H,K,K1,K2,C,BE,S,thickness',
                        help='comma separated curvature function list to process')
    parser.add_argument('--statFunc', '-f',
                        dest='statFunc',
                        action='store',
                        default='ptile-raw',
                        help='comma separated statistical function list to process')
    parser.add_argument('--hemi', '-m',
                        dest='hemi',
                        action='store',
                        default='lh,rh',
                        help='comma separated hemisphere list to process')
    parser.add_argument('--schedulerStdOutDir',
                        dest='schedulerStdOutDir',
                        action='store',
                        default='~/scratch',
                        help='top level directory containing stdout from scheduled jobs')
    parser.add_argument('--schedulerStdErrDir',
                        dest='schedulerStdErrDir',
                        action='store',
                        default='~/scratch',
                        help='top level directory containing stderr from scheduled jobs')

    args = parser.parse_args()

    # A generic "shell"
    OSshell = crun.crun()
    OSshell.echo(False)
    OSshell.echoStdOut(args.subShellFlush)
    # OSshell.runCmd(not args.subShellDontRunCmd)
    OSshell.detach(False)

    roicontroller = FNNDSC_ROICONTROLLER(
                        subjectList             = args.l_subj,
                        sampleList              = args.sample,
                        annot                   = args.annotation,
                        outDir                  = args.outDir,
                        workingDir              = args.workingDir,
                        stages                  = args.stages,
                        pval                    = args.pval,
                        group                   = args.group,
                        surface                 = args.surface,
                        statFunc                = args.statFunc,
                        curvFunc                = args.curvFunc,
                        hemi                    = args.hemi,
                        sampleThreshold         = args.threshold,
                        schedulerStdOutDir      = args.schedulerStdOutDir,
                        schedulerStdErrDir      = args.schedulerStdErrDir,
                        logTo                   = 'ROIcontroller.log',
                        syslog                  = True,
                        logTee                  = True)

    roicontroller.verbosity(args.verbosity)
    pipeline    = roicontroller.pipeline()
    pipeline.poststdout(True)
    pipeline.poststderr(True)
    roicontroller.topDir(os.getcwd())

    stage0 = stage.Stage(
                        name            = 'ROIcontroller-0-init',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIcontroller-0-init.log',
                        logTee          = True,
                        )
    def f_stage0callback(**kwargs):
        '''
        STAGE 0
        -------

        Build output directory trees for the analysis.

        '''
        str_cwd         =  os.getcwd()
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        log                    = stage.log()
        b_alreadyResampled     = False
        b_outAnnotDirExist     = False

        for pipeline._str_annotation in pipeline.l_annotation():
            for pipeline._str_subsample in pipeline.l_subsample():
                b_parentContainsSamples = False
                str_outDir    = "%s%s/groupCurvAnalysis/%s" % \
                    (   pipeline.outDir(),
                        pipeline._str_subsample,
                        pipeline._str_annotation    )
                str_outParentDir = "/".join(str_outDir.split('/')[0:-2])
                log('Building output tree for "%s", annotation "%s"\n' % \
                    (   str_outParentDir,
                        pipeline._str_annotation))

                b_outAnnotDirExist = os.path.isdir(str_outDir)
                if not b_outAnnotDirExist:
                    misc.mkdir(str_outDir)

                # Check if outParentDir contains samples
                l_numer = l_filterNumeric(os.listdir(str_outParentDir))
                if len(l_numer): b_parentContainsSamples = True
        stage.exitCode(0)
        return True
    stage0.def_stage(f_stage0callback, obj=stage0, pipe=roicontroller)

    stage1 = stage.Stage(
                        name            = 'ROIcontroller-1-resample',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIcontroller-1-resample.log',
                        logTee          = True
                        )
    stage1.def_preconditions(stage0._f_callCount)
    def f_stage1callback(**kwargs):
        '''
        STAGE 1
        -------

        Run the resampler (intelligently) across all the sample experiment
        directories.

        '''
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        os.chdir(pipeline.topDir())
        log                    = stage.log()

        for pipeline._str_annotation in pipeline.l_annotation():
            for pipeline._str_subsample in pipeline.l_subsample():
                b_parentContainsSamples = False
                str_outDir    = "%s%s" % \
                    (
                        pipeline.outDir(),
                        pipeline._str_subsample,
                    )
                log('outDir = %s\n' % str_outDir)

                # Check if outParentDir contains samples
                l_numer = l_filterNumeric(os.listdir(str_outDir))
                log('numerical directory count = %d\n' % len(l_numer))
                if len(l_numer):
                    b_parentContainsSamples = True
                    str_stages = "12"
                else:
                    str_stages = "012"
                str_subjList = ' '.join(pipeline._l_subj)
                str_execCmd = "%s -o %s -r %s -s %s -a %s %s" % \
                    (
                        "/neuro/users/rudolphpienaar/src/devel/roi_tag/rsample.py",
                        str_outDir,
                        pipeline._str_threshold,
                        str_stages,
                        pipeline._str_annotation,
                        str_subjList
                    )
                log("Shell command = %s\n" % str_execCmd)
                OSshell(str_execCmd)
        stage.exitCode(0)
        return True
    stage1.def_stage(f_stage1callback, obj=stage1, pipe=roicontroller)

    stage2 = stage.Stage(
                        name            = 'ROIcontroller-2-curvatureAnalysis',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIcontroller-2-curvatureAnalysis.log',
                        logTee          = True
                        )
    stage2.def_preconditions(stage0._f_callCount)
    def f_stage2callback(**kwargs):
        '''
        STAGE 2
        -------

        Run the (new) 2D cloud-based curvature analysis.

        '''
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        os.chdir(pipeline.topDir())
        log                    = stage.log()
        log('Stage 2 -- running from dir %s\n' % pipeline.topDir())

        shellCluster = crun.crun_hpc_mosix(remoteUser='rudolphpienaar',
                                           remoteHost='rc-majesty',
                                           schedulerStdOutDir='%s-%s' % (pipeline.schedulerStdOutDir(), '-stage-2'),
                                           schedulerStdErrDir='%s-%s' % (pipeline.schedulerStdErrDir(), '-stage-2'))
        for pipeline._str_annotation in pipeline.l_annotation():
            for pipeline._str_subsample in pipeline.l_subsample():
                '''
                Read the ROI.lst file containing the various annotation ROIs,
                and for each element, start a separate deviation analysis
                process.
                '''

                str_outDir    = "%s%s/groupCurvAnalysis/%s" % \
                    (   pipeline.outDir(),
                        pipeline._str_subsample,
                        pipeline._str_annotation    )
                os.chdir("%s/%s" % (pipeline.topDir(), str_outDir))

                with open('ROI.lst') as f: s = f.read()
                l_ROI = s.strip().split(' ')
                for str_ROI in l_ROI:
                    for pipeline._str_hemi in pipeline.l_hemisphere():
                        for pipeline._str_surface in pipeline.l_surface():
                            str_execCmd = '%s/../../sh/%s -D %s/%s -h %s -s %s %s' % \
                            (
                                pipeline.topDir(),
                                './ptile_deviation_analysis-raw.bash',
                                pipeline.topDir(), str_outDir,
                                pipeline._str_hemi,
                                pipeline._str_surface,
                                str_ROI
                            )
                            print(os.getcwd())
                            log("Shell command = %s\n" % str_execCmd)
                            #OSshell.echo(True)
                            #OSshell.echoStdOut(True)
                            #OSshell.echoStdErr(True)
                            #OSshell.waitForChild(True)
                            #OSshell(str_execCmd)
                            shellCluster.waitForChild(True)
                            shellCluster(str_execCmd)
                            # Sleep for a second so as not to overload the head node
                            # in tight scheduling loops
                            time.sleep(1)
        shellCluster.blockOnChild()
        stage.exitCode(0)
        return True
    stage2.def_stage(f_stage2callback, obj=stage2, pipe=roicontroller)


    stage3 = stage.Stage(
                        name            = 'ROIcontroller-3-dty_analysis',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIcontroller-3-dty_analysis.log',
                        logTee          = True
                        )
    stage3.def_preconditions(stage0._f_callCount)
    def f_stage3callback(**kwargs):
        '''
        STAGE 3
        -------

        Run the overlap/density analysis in the cloud space.

        '''
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        os.chdir(pipeline.topDir())
        log                    = stage.log()
        log('Stage 3 -- running from dir %s\n' % pipeline.topDir())
        shellCluster = crun.crun_hpc_mosix(remoteUser='rudolphpienaar',
                                           remoteHost='rc-majesty',
                                           schedulerStdOutDir='%s-%s' % (pipeline.schedulerStdOutDir(), '-stage-3'),
                                           schedulerStdErrDir='%s-%s' % (pipeline.schedulerStdErrDir(), '-stage-3'))

        for pipeline._str_annotation in pipeline.l_annotation():
            for pipeline._str_subsample in pipeline.l_subsample():
                for pipeline._str_pval in pipeline.l_pval():
                    str_outDir    = "%s%s/groupCurvAnalysis/%s" % \
                        (   pipeline.outDir(),
                            pipeline._str_subsample,
                            pipeline._str_annotation    )
                    os.chdir("%s/%s" % (pipeline.topDir(), str_outDir))
                    str_groups  = ','.join(pipeline.l_group())
                    str_curv    = ','.join(pipeline.l_curvFunc())
                    for pipeline._str_hemi in pipeline.l_hemisphere():
                        for pipeline._str_surface in pipeline.l_surface():
                            str_execCmd = '%s/../../sh/%s -D %s/%s -p %s -g %s -h %s -c %s -s %s' % \
                                        (
                                            pipeline.topDir(),
                                            "./dty_analyze-here.bash",
                                            pipeline.topDir(), str_outDir,
                                            pipeline._str_pval,
                                            str_groups,
                                            pipeline._str_hemi,
                                            str_curv,
                                            pipeline._str_surface
                                        )
                            log("Current dir = %s\n" % os.getcwd())
                            log("Shell command = %s\n" % str_execCmd)
                            shellCluster.waitForChild(True)
                            shellCluster(str_execCmd)
        shellCluster.blockOnChild()
        stage.exitCode(0)
        return True
    stage3.def_stage(f_stage3callback, obj=stage3, pipe=roicontroller)


    stage4 = stage.Stage(
                        name            = 'ROIcontroller-4-roi_tag',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIcontroller-4-roi_tag.log',
                        logTee          = True
                        )
    stage4.def_preconditions(stage0._f_callCount)
    def f_stage4callback(**kwargs):
        '''
        STAGE 4
        -------

        Run the ROI tagging.

        '''
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        os.chdir(pipeline.topDir())
        log                    = stage.log()
        log('Stage 4 -- running from dir %s\n' % pipeline.topDir())
        shellCluster = crun.crun_hpc_mosix(remoteUser='rudolphpienaar',
                                           remoteHost='rc-majesty',
                                           schedulerStdOutDir='%s-%s' % (pipeline.schedulerStdOutDir(), '-stage-4'),
                                           schedulerStdErrDir='%s-%s' % (pipeline.schedulerStdErrDir(), '-stage-4'))

        for pipeline._str_annotation in pipeline.l_annotation():
            for pipeline._str_subsample in pipeline.l_subsample():
                str_outDir    = "%s%s/groupCurvAnalysis/%s" % \
                    (   pipeline.outDir(),
                        pipeline._str_subsample,
                        pipeline._str_annotation    )
                os.chdir("%s/%s" % (pipeline.topDir(), str_outDir))
                str_pval    = ','.join(pipeline.l_pval())
                str_groups  = ','.join(pipeline.l_group())
                str_surface = ','.join(pipeline.l_surface())
                str_statFunc= ','.join(pipeline.l_statFunc())
                str_hemi    = ','.join(pipeline.l_hemisphere())
                str_curv    = ','.join(pipeline.l_curvFunc())

                str_execCmd = '%s --clobber --workingDir %s/%s --curvFunc %s --pval %s --group %s --surface %s --statFunc %s --hemi %s --stages 01234 --schedulerStdOutDir=%s --schedulerStdErrDir=%s $(cat ./ROI.lst)' % \
                    (
#                        pipeline.topDir(), str_outDir,
                        "/neuro/users/rudolphpienaar/src/devel/roi_tag/roi_tag.py",
                        pipeline.topDir(), str_outDir,
                        str_curv,
                        str_pval,
                        str_groups,
                        str_surface,
                        str_statFunc,
                        str_hemi,
                        '/neuro/users/rudolphpienaar/scratch/%s-roi_tag.jobout' % os.getpid(),
                        '/neuro/users/rudolphpienaar/scratch/%s-roi_tag.joberr' % os.getpid()
                    )
                log("Shell command = %s\n" % str_execCmd)
                shellCluster.waitForChild(True)
                shellCluster(str_execCmd)
        shellCluster.blockOnChild()
        stage.exitCode(0)
        return True
    stage4.def_stage(f_stage4callback, obj=stage4, pipe=roicontroller)


    roicontrollerlog = roicontroller.log()
    roicontrollerlog('INIT: %s\n' % ' '.join(sys.argv))
    roicontroller.stage_add(stage0)
    roicontroller.stage_add(stage1)
    roicontroller.stage_add(stage2)
    roicontroller.stage_add(stage3)
    roicontroller.stage_add(stage4)
    roicontroller.initialize()


    roicontroller.run()
