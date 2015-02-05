#!/usr/bin/env python

'''

    'roi_tag.py' analyzes the output directories of a curvature analysis run
    and groups "tagged" ROI label files in a manner that allows for
    simplified display/projection on brain surfaces.

'''


import  os
import  sys
import  argparse
import  tempfile, shutil
from    _common import systemMisc       as misc
from    _common import crun

import  error
import  message
import  stage
import  fnndsc  as base

class FNNDSC_ROITAG(base.FNNDSC):
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
        'workingDirNotExist': {
            'action'        : 'attempting to access the <workingDir>, ',
            'error'         : 'a system error was encountered. Does the directory exist?',
            'exitCode'      : 16},
    }

    def l_pval(self):
        return self._l_pval

    def l_roi(self):
        return self._l_ROI

    def l_hemisphere(self):
        return self._l_hemi

    def l_surface(self):
        return self._l_surface

    def l_statFunc(self):
        return self._l_statFunc

    def l_group(self):
        return self._l_group

    def l_curvFunc(self):
        return self._l_curvFunc

    def pval(self):
        return self._str_pval

    def topDir(self, *args):
        if len(args):
            self._topDir    = args[0]
        else:
            return self._topDir

    def dirSpec(self):
        '''
        Return the dirSpec based on internal pipeline._str_* variables
        '''
        return '%s/%s/%s/%s/%s/%s' % (  self.outDir(),
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
        return '%s/%s/%s/%s' % (    self._str_pval,
                                    self._str_statFunc,
                                    self._str_surface,
                                    self._str_hemi)

    def namespec(self, *args):
        '''
        Return the namespec based on internal pipeline._str_* variables.
        '''
        str_sep = "-"
        if len(args): str_sep = args[0]
        return '%s%s%s%s%s%s%s%s%s' % ( self._str_group,     str_sep,
                                        self._str_pval,      str_sep,
                                        self._str_statFunc,  str_sep,
                                        self._str_surface,   str_sep,
                                        self._str_hemi)

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

    def roi(self):
        return self._str_roi

    def surface(self):
        return self._str_surface

    def hemi(self):
        return self._str_hemi

    def statFunc(self):
        return self._str_statFunc

    def curvFunc(self):
        return self._str_curvFunc

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

    def clobber(self, *args):
        if len(args):
            self._b_clobber = args[0]
        else:
            return self._b_clobber

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

        self._l_pval                    = []
        self._l_group                   = []
        self._l_surface                 = []
        self._l_statFunc                = []
        self._l_curvFunc                = []
        self._l_hemi                    = []

        self._outDir                    = ''
        self._workingDir                = ''
        self._stageslist                = '12'

        # Internal tracking vars
        self._str_pval                  = ''
        self._str_group                 = ''
        self._str_roi                   = ''
        self._str_hemi                  = ''
        self._str_surface               = ''
        self._str_statFunc              = ''
        self._str_curvFunc              = ''

        self._topDir                    = ''

        # Scheduler std out/err dirs
        self._str_schedulerStdOutDir    = '~/scratch'
        self._str_schedulerStdErrDir    = '~/scratch'

        self._b_clobber                 = False

        for key, value in kwargs.iteritems():
            if key == 'ROIList':        self._l_ROI             = value
            if key == 'outDir':         self._outDir            = value
            if key == 'workingDir':     self._workingDir        = value
            if key == 'stages':         self._stageslist        = value
            if key == 'pval':           self._l_pval            = value.split(',')
            if key == 'group':          self._l_group           = value.split(',')
            if key == 'surface':        self._l_surface         = value.split(',')
            if key == 'statFunc':       self._l_statFunc        = value.split(',')
            if key == 'curvFunc':       self._l_curvFunc        = value.split(',')
            if key == 'hemi':           self._l_hemi            = value.split(',')
            if key == 'schedulerStdOutDir':     self._str_schedulerStdOutDir = value
            if key == 'schedulerStdErrDir':     self._str_schedulerStdErrDir = value

        if not os.path.isdir(self._workingDir): errorFatal(self, 'workingDirNotExist')


    def initialize(self):
        '''
        This method provides some "post-constructor" initialization. It is
        typically called after the constructor and after other class flags
        have been set (or reset).

        '''

        # Set the stages
        self._pipeline.stages_canRun(False)
        lst_stages = list(self._stageslist)
        for index in lst_stages:
            stage = self._pipeline.stage_get(int(index))
            stage.canRun(True)

    def run(self):
        '''
        The main 'engine' of the class.

        '''
        base.FNNDSC.run(self)


def synopsis(ab_shortOnly = False):
    scriptName = os.path.basename(sys.argv[0])
    shortSynopsis =  '''
    SYNOPSIS

            %s                                            \\
                            [--stages <stages>]             \\
                            [-o|--outDir <outputRootDir>]   \\
                            [-w|--workingDir <workingDir>]  \\
                            [-v|--verbosity <verboseLevel>] \\
                            [-s|--stages <stages>]          \\
                            [-p|--pval <pvalCutoffList>]    \\
                            [-g|--group <groupList>]        \\
                            [-S|--surface <surfaceList>]    \\
                            [-f|--statFunc <statFuncList>]  \\
                            [-c|--curvFunc <curvFuncList>]  \\
                            [-m|--hemi <hemiList>]          \\
                            [--schedulerStdOutDir <dir>]    \\
                            [--schedulerStdErrDir <dir>]    \\
                            <ROI_1> <ROI_2> ... <ROI_N>
    ''' % scriptName

    description =  '''
    DESCRIPTION

        `%s' creates ROI label summaries. Once a curvature and statistical
        geometric density analysis has been performed, this script will create
        a list of surface label files that have been "tagged" by the analysis
        as having significant group separation.

        Results are organized per group, and further categorized into ROIs that
        are

            * pval only significant
            * statistically geometric only significant
            * pval *and* statistical geometric significant

        thickness values are a special case and saved to their own tree.

        This tree structure is thus:

        <group>
          +
          +--<pval>
            +
            +--<statFunc>
              +
              +--<surface>
                +
                +--<hemi>
                  +
                  |
                  +-- K
                  | +-- pval
                  | +-- statFunc
                  | +-- pval-statFunc
                  +-- K1
                  | +-- pval
                  | +-- statFunc
                  | +-- pval-statFunc
                  +-- K2
                  | +-- pval
                  | +-- statFunc
                  | +-- pval-statFunc

                        ...

                  +-- thickness
                    +-- pval
                    +-- statFunc
                    +-- pval-statFunc

        Each of the above subdirectories will contain a list of label ROI files
        that satisfy the named constraint for the given curvature or thickness
        distribution.

    ARGS

       --stages <stages>
       The stages to execute. This is specified in a string, such as '1234'
       which would imply stages 1, 2, 3, and 4.

       The special keyword 'all' can be used to turn on all stages.

       --pval <pvalCutoffList>
       The pval cutoffs to consider. In practice, this is always 'le1,le5'

       --group <groupList>
       The group list to process.

       --surface <surfaceList>
       The surface list to process. In practice, 'smoothwm,pial'.

       --statFunc <statFuncList>
       The statistical functional data to analyze. Typically
       'ptile-raw,ptile-convex'.

       --curvFunc <curvFuncList>
       The curvature functions to analyze. Typically
       'K,K1,K2,H,BE,S,C,thickness'.

       --hemi <hemiList>
       The hemispheres to process. In practice, this is always 'lh,rh'.

       --workingDir <workingDir>
       The working directory for the script.

       --output <outputRootDir>
       The top level directory name to contain results. The fully qualified
       output dir is <workingDir>/<outputDir>

       --clobber
       A boolean flag. If true, will not delete existing output directories
       but simply add more results down existing trees, clobbering existing
       files if they already exist. If not specified then existing output
       trees are deleted. This assures that a given run contains only data
       from that run.

       Note that running the same experiment multiple times with "--clobber"
       will *grow* resultant ROI label files!

       For a distributed experiment, delete *all* existing ROItag trees
       *before* running the experiment!

       <ROI_1> <ROI_2> ... <ROI_N>
       The ROIs to process.

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

    parser.add_argument('l_roi',
                        metavar='ROIs', nargs='+',
                        help='ROIs to process')
    parser.add_argument('--verbosity', '-v',
                        dest='verbosity',
                        action='store',
                        default=0,
                        help='verbosity level')
    parser.add_argument('--output', '-o',
                        dest='outDir',
                        action='store',
                        default='ROItag',
                        help='output root directory')
    parser.add_argument('--workingDir', '-w',
                        dest='workingDir',
                        action='store',
                        default='./',
                        help='output working directory')
    parser.add_argument('--clobber',
                        dest='clobber',
                        action='store_true',
                        default=False,
                        help='if specified, do not erase existing output dir if found.')
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
    parser.add_argument('--statFunc', '-f',
                        dest='statFunc',
                        action='store',
                        default='ptile-raw',
                        help='comma separated statistical function list to process')
    parser.add_argument('--curvFunc', '-c',
                        dest='curvFunc',
                        action='store',
                        default='H,K,K1,K2,C,BE,S,thickness',
                        help='comma separated curvature function list to process')
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
    OSshell.echoStdOut(False)
    OSshell.detach(False)
    OSshell.waitForChild(True)

    roitag = FNNDSC_ROITAG(
                        ROIList                 = args.l_roi,
                        outDir                  = '%s-%s' % (args.outDir, args.curvFunc),
                        workingDir              = args.workingDir,
                        stages                  = args.stages,
                        pval                    = args.pval,
                        group                   = args.group,
                        surface                 = args.surface,
                        curvFunc                = args.curvFunc,
                        statFunc                = args.statFunc,
                        hemi                    = args.hemi,
                        schedulerStdOutDir      = args.schedulerStdOutDir,
                        schedulerStdErrDir      = args.schedulerStdErrDir,
                        logTo                   = '%s/ROItag.log' % args.workingDir,
                        syslog                  = True,
                        logTee                  = True)

    roitag.clobber(args.clobber)
    roitag.verbosity(args.verbosity)
    pipeline    = roitag.pipeline()
    pipeline.poststdout(True)
    pipeline.poststderr(True)
    roitag.topDir(os.getcwd())

    stage0 = stage.Stage(
                        name            = 'ROItag-0-init',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/ROItag-0-init.log' % args.workingDir,
                        logTee          = True,
                        )
    def f_stage0callback(**kwargs):
        str_cwd         =  os.getcwd()
        for key, val in kwargs.iteritems():
            if key == 'roi':    l_roi       = val
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        str_cwd         = pipeline.workingDir()
        log             = stage.log()
        lst_pval        = pipeline.l_pval()
        lst_group       = pipeline.l_group()
        lst_hemi        = pipeline.l_hemisphere()
        lst_surface     = pipeline.l_surface()
        lst_statFunc    = pipeline.l_statFunc()
        lst_roi         = pipeline.l_roi()
        lst_ctype       = pipeline.l_curvFunc()
        lst_ctype       = ['thickness', 'curv']
        lst_leaf        = ['pval', 'statFunc', 'pval-only', 'statFunc-only', 'pval_N_statFunc']
        os.chdir(str_cwd)
        if os.path.isdir(pipeline.outDir()) and not pipeline.clobber():
            log('Existing outDir tree found... deleting...\n')
            shutil.rmtree(pipeline.outDir())
        OSshell('mkdir -p %s' % pipeline.outDir())
        os.chdir(pipeline.outDir())
        pipeline.outDir(os.getcwd())
        os.chdir(str_cwd)
        if OSshell.exitCode() != 0: error.fatal(pipeline, 'outDirNotCreate')

        for pipeline._str_group in lst_group:
            for pipeline._str_pval in lst_pval:
                for pipeline._str_statFunc in lst_statFunc:
                    for pipeline._str_surface in lst_surface:
                        for pipeline._str_hemi in lst_hemi:
                            for ctype in lst_ctype:
                                for leaf in lst_leaf:
                                    OSshell('mkdir -p %s/%s/%s' % (
                                                            pipeline.dirSpec(),
                                                            ctype,
                                                            leaf))
        stage.exitCode(0)
        return True

    stage0.def_stage(f_stage0callback, roi=args.l_roi, obj=stage0, pipe=roitag)
#    stage0.def_postconditions(f_stageShellExitCode, obj=stage0)

    stage1 = stage.Stage(
                        name            = 'ROItag-1-sort-pval',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/ROItag-1-sort-pval.log'  % args.workingDir,
                        logTee          = True
                        )
    # stage1.def_preconditions(stage0.def_postconditions()[0], **stage0.def_postconditions()[1])
    stage1.def_preconditions(lambda **x: True)
    def f_stage1callback(**kwargs):
        ''' Builds list of p-val tagged hits '''
        #str_cwd         =  os.getcwd()
        l_roi           = []

        for key, val in kwargs.iteritems():
            if key == 'roi':    l_roi       = val
            if key == 'pipe':   pipeline    = val
            if key == 'obj':    stage       = val
        str_cwd         = pipeline.workingDir()
        log             = stage.log()
        lst_statFunc    = pipeline.l_statFunc()
        lst_group       = pipeline.l_group()
        lst_pval        = pipeline.l_pval()
        for roi in l_roi:
            log('Processing %s\n' % roi)
            for pipeline._str_statFunc in lst_statFunc:
                for pipeline._str_pval in lst_pval:
                    for pipeline._str_surface in pipeline.l_surface():
                        for pipeline._str_hemi in pipeline.l_hemisphere():
                            str_summaryFileDir = '%s/%s/dty_analyze-1-summaryTables-%s-%s' % \
                                                    (roi, pipeline._str_statFunc, pipeline._str_pval,
                                                     args.curvFunc)
                            # print(str_summaryFileDir)
                            os.chdir('%s/%s' % (pipeline.workingDir(), str_summaryFileDir))
                            for pipeline._str_group in lst_group:
                                str_tagFile = 'p-%s-%s-%s-%s' % \
                                    (pipeline._str_group, pipeline._str_pval,
                                     pipeline._str_hemi,  pipeline._str_surface)
                                log('Current dir is %s\n' % os.getcwd())
                                if os.path.getsize(str_tagFile):
                                    for hit in [line.strip() for line in open(str_tagFile, 'r')]:
                                        l_fileSpec  = os.path.basename(hit).split('-')
                                        pipeline._str_hemi      = l_fileSpec[2]
                                        pipeline._str_func      = l_fileSpec[3]
                                        pipeline._str_surface   = l_fileSpec[-4]
                                        if pipeline._str_func != 'thickness':
                                            ctype = 'curv'
                                        else:
                                            ctype = 'thickness'
                                        # print(os.getcwd())
                                        # print(l_fileSpec)
                                        # log('hit = %s\n' % hit)

                                        curvMatch = [s for s in pipeline.l_curvFunc() if pipeline._str_func == s]
                                        if len(curvMatch):

                                            log('Logging hit for group %s, %s.%s.%s, p-val %s, ctype %s\n' % \
                                                (   pipeline._str_group,
                                                    pipeline._str_hemi,
                                                    pipeline._str_func,
                                                    pipeline._str_surface,
                                                    pipeline._str_pval,
                                                    ctype))
                                            try:
                                                misc.file_writeOnce('%s/%s/pval/%s' % \
                                                            (pipeline.dirSpec(), ctype, roi), '%s\n' % hit, mode='a')
                                            except:
                                                log('Ignoring %s-%s\n' % (pipeline._str_hemi, pipeline._str_surface))

                            os.chdir(str_cwd)
        return True
    stage1.def_stage(f_stage1callback, roi=args.l_roi, obj=stage1, pipe=roitag)
    # stage1.def_postconditions(f_blockOnScheduledJobs, obj=stage1,
    #                           blockProcess    = 'mris_label_calc')


    stage2 = stage.Stage(
                        name            = 'ROItag-2-sort-areaDensity',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/ROItag-2-sort-areaDensity.log'  % args.workingDir,
                        logTee          = True
                        )
    # stage2.def_preconditions(stage1.def_postconditions()[0], **stage1.def_postconditions()[1])
    stage2.def_preconditions(lambda **x: True)
    def f_stage2callback(**kwargs):
        ''' Builds list of statFunc tagged hits '''
        str_cwd         =  os.getcwd()
        l_roi           = []

        for key, val in kwargs.iteritems():
            if key == 'roi':    l_roi       = val
            if key == 'pipe':   pipeline    = val
            if key == 'obj':    stage       = val
        str_cwd         = pipeline.workingDir()
        log             = stage.log()
        lst_statFunc    = pipeline.l_statFunc()
        lst_group       = pipeline.l_group()
        lst_pval        = pipeline.l_pval()
        for roi in l_roi:
            log('Processing %s\n' % roi)
            for pipeline._str_statFunc in lst_statFunc:
                for pipeline._str_pval in lst_pval:
                    for pipeline._str_surface in pipeline.l_surface():
                        for pipeline._str_hemi in pipeline.l_hemisphere():
                            str_summaryFileDir = '%s/%s/dty_analyze-3-cutoffs-%s-%s' % \
                                                    (roi, pipeline._str_statFunc, pipeline._str_pval,
                                                     args.curvFunc)
                            # print(str_summaryFileDir)
                            os.chdir(str_summaryFileDir)
                            for pipeline._str_group in lst_group:
                                str_tagFile = 'separate-%s-AreaDensity-%s-%s.txt' % \
                                    (pipeline._str_group, pipeline._str_hemi, pipeline._str_surface)
                                if os.path.isfile(str_tagFile):
                                    for hit in [line.strip() for line in open(str_tagFile, 'r')]:
                                        hitFile = hit.split()[2]
                                        l_fileSpec  = os.path.basename(hitFile).split('-')
                                        pipeline._str_hemi      = l_fileSpec[2]
                                        pipeline._str_func      = l_fileSpec[3]
                                        pipeline._str_surface   = l_fileSpec[-3]
                                        if pipeline._str_func != 'thickness':
                                            ctype = 'curv'
                                        else:
                                            ctype = 'thickness'
                                        # print(os.getcwd())
                                        # print(l_fileSpec)
                                        # log('hit = %s\n' % hit)

                                        curvMatch = [s for s in pipeline.l_curvFunc() if pipeline._str_func == s]
                                        if len(curvMatch):

                                            log('Logging hit for group %s, %s.%s.%s, p-val theshold %s, ctype %s\n' % \
                                                (   pipeline._str_group,
                                                    pipeline._str_hemi,
                                                    pipeline._str_func,
                                                    pipeline._str_surface,
                                                    pipeline._str_pval,
                                                    ctype))
                                            try:
                                                misc.file_writeOnce('%s/%s/statFunc/%s' % \
                                                            (pipeline.dirSpec(), ctype, roi), '%s\n' % hit, mode='a')
                                            except:
                                                log('Ignoring %s-%s\n' % (pipeline._str_hemi, pipeline._str_surface))

                            os.chdir(str_cwd)
        return True
    stage2.def_stage(f_stage2callback, roi=args.l_roi, obj=stage2, pipe=roitag)
    # stage2.def_postconditions(f_blockOnScheduledJobs, obj=stage2,
    #                           blockProcess    = 'mris_label_calc')


    def tcl_append(str_prefix, str_suffix):
        '''
        Append some text to each tcl file to save snapshots of different brain aspects.
        '''
        str_content = '''
        # Initial medial view
        read_binary_curv
        redraw
        save_tiff %s-medial-%s.tiff

        # inferior view
        rotate_brain_x 90
        redraw
        save_tiff %s-inferior-%s.tiff

        # superior view
        rotate_brain_x -180
        redraw
        save_tiff %s-superior-%s.tiff

        # reset
        rotate_brain_x 90

        # frontal view
        rotate_brain_y 90
        redraw
        save_tiff %s-frontal-%s.tiff

        # lateral view
        rotate_brain_y 90
        redraw
        save_tiff %s-lateral-%s.tiff

        # dorsal view
        rotate_brain_y 90
        redraw
        save_tiff %s-dorsal-%s.tiff

        exit 0
        ''' % ( str_prefix, str_suffix,
                str_prefix, str_suffix,
                str_prefix, str_suffix,
                str_prefix, str_suffix,
                str_prefix, str_suffix,
                str_prefix, str_suffix)
        return str_content

    stage3 = stage.Stage(
                        name            = 'ROItag-3-sort-pval-statFunc',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/ROItag-3-sort-pval-statFunc.log'  % args.workingDir,
                        logTee          = True
                        )
    # stage3.def_preconditions(stage1.def_postconditions()[0], **stage1.def_postconditions()[1])
    stage3.def_preconditions(lambda **x: True)
    def f_stage3callback(**kwargs):
        ''' Sorts intersect/xor of classes '''
        str_cwd         =  os.getcwd()
        for key, val in kwargs.iteritems():
            if key == 'roi':    l_roi       = val
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        str_cwd         = pipeline.workingDir()
        log             = stage.log()
        lst_pval        = pipeline.l_pval()
        lst_group       = pipeline.l_group()
        lst_hemi        = pipeline.l_hemisphere()
        lst_surface     = pipeline.l_surface()
        lst_statFunc    = pipeline.l_statFunc()
        lst_roi         = pipeline.l_roi()
        lst_ctype       = ['thickness', 'curv']
        lst_leaf        = ['pval', 'statFunc', 'pval-only', 'statFunc-only', 'pval_N_statFunc']

        for pipeline._str_group in lst_group:
            for pipeline._str_pval in lst_pval:
                for pipeline._str_statFunc in lst_statFunc:
                    for pipeline._str_surface in lst_surface:
                        for pipeline._str_hemi in lst_hemi:
                            for ctype in lst_ctype:
                                log('Processing %s\n' % pipeline.dirSpecPartial())
                                os.chdir('%s/%s/%s' % (pipeline.dirSpec(), ctype, 'pval'))
                                l_pval      = OSshell('ls')[0].strip().split('\n')
                                os.chdir('%s/%s/%s' % (pipeline.dirSpec(), ctype, 'statFunc'))
                                l_statFunc  = OSshell('ls')[0].strip().split('\n')
                                os.chdir('../')
                                if l_pval == ['']:      sp = set()
                                else:                   sp = set(l_pval)
                                if l_statFunc == ['']:  sf = set()
                                else:                   sf = set(l_statFunc)
                                U   = sp.union(sf)
                                I   = sp.intersection(sf)
                                pO  = sp - I
                                fO  = sf - I
                                for roi in list(I):
                                    log('Intersection ROI: %s\n' % roi)
                                    shutil.copyfile('statFunc/%s' % roi, 'pval_N_statFunc/%s' % roi )
                                for roi in list(pO):
                                    log('p-val only ROI: %s\n' % roi)
                                    shutil.copyfile('pval/%s' % roi, 'pval-only/%s' % roi)
                                for roi in list(fO):
                                    log('statFunc only ROI: %s\n' % roi)
                                    shutil.copyfile('statFunc/%s' % roi, 'statFunc-only/%s' % roi)
        stage.exitCode(0)
        return True
    stage3.def_stage(f_stage3callback, roi=args.l_roi, obj=stage3, pipe=roitag)
    # stage2.def_postconditions(f_blockOnScheduledJobs, obj=stage2,
    #                           blockProcess    = 'mris_label_calc')

    stage4 = stage.Stage(
                        name            = 'ROItag-4-labelReader',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/ROItag-4-labelReader.log'  % args.workingDir,
                        logTee          = True
                        )
    # stage4.def_preconditions(stage1.def_postconditions()[0], **stage1.def_postconditions()[1])
    stage4.def_preconditions(lambda **x: True)
    def f_stage4callback(**kwargs):
        str_cwd         =  os.getcwd()
        for key, val in kwargs.iteritems():
            if key == 'roi':    l_roi       = val
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        str_cwd         = pipeline.workingDir()
        log             = stage.log()
        lst_pval        = pipeline.l_pval()
        lst_group       = pipeline.l_group()
        lst_hemi        = pipeline.l_hemisphere()
        lst_surface     = pipeline.l_surface()
        lst_statFunc    = pipeline.l_statFunc()
        lst_roi         = pipeline.l_roi()
        lst_ctype       = ['thickness', 'curv']
        lst_leaf        = ['pval', 'statFunc', 'pval-statFunc']
        shellCluster = crun.crun_hpc_mosix(remoteUser='rudolphpienaar',
                                           remoteHost='rc-majesty',
                                           schedulerStdOutDir='%s-%s' % (pipeline.schedulerStdOutDir(), '-stage-4'),
                                           schedulerStdErrDir='%s-%s' % (pipeline.schedulerStdErrDir(), '-stage-4'))
        str_curvFunc = ','.join(pipeline.l_curvFunc())
        log('Group list: %s\n' % lst_group)
        for pipeline._str_group in lst_group:
            for pipeline._str_pval in lst_pval:
                for pipeline._str_statFunc in lst_statFunc:
                    for pipeline._str_surface in lst_surface:
                        for pipeline._str_hemi in lst_hemi:
                            for ctype in lst_ctype:
                                if ctype == 'thickness': str_curvFunc = 'thickness'
                                else:
                                    l_curv = list(pipeline.l_curvFunc())
                                    log('curvFunc = %s\n' % l_curv)
                                    if any('thickness' in s for s in l_curv):
                                        l_curv.remove('thickness')
                                    str_curvFunc = ','.join(l_curv)
                                log('Processing %s-%s-%s\n' % \
                                    (pipeline._str_group, pipeline.dirSpecPartial(), ctype))
                                os.chdir('%s/%s/%s' % (pipeline.dirSpec(), ctype, 'pval-only'))
                                sortCmd     = "/bin/ls -l | sort -n -k 5 | awk '{print $9}'"
                                l_pval      = OSshell(sortCmd)[0].strip().split('\n')
                                l_pval      = filter(None, l_pval)
                                os.chdir('%s/%s/%s' % (pipeline.dirSpec(), ctype, 'statFunc-only'))
                                l_statFunc  = OSshell(sortCmd)[0].strip().split('\n')
                                l_statFunc  = filter(None, l_statFunc)
                                os.chdir('%s/%s/%s' % (pipeline.dirSpec(), ctype, 'pval_N_statFunc'))
                                l_pNs       = OSshell(sortCmd)[0].strip().split('\n')
                                l_pNs       = filter(None, l_pNs)
                                os.chdir('../')
                                labelCount      = 0
                                rLabelCount     = 0
                                colorMult       = 9
                                misc.file_writeOnce('pval-only.freeview.txt', '#!/bin/bash\n\nSUBJNAME=$1; SURF=$2\nexport SUBJECT=$SUBJECTS_DIR/$SUBJNAME\nfreeview --surface $SUBJECT/surf/%s.$SURF' % pipeline._str_hemi)
                                misc.file_writeOnce('statFunc-only.freeview.txt', '#!/bin/bash\n\nSUBJNAME=$1; SURF=$2\nexport SUBJECT=$SUBJECTS_DIR/$SUBJNAME\nfreeview --surface $SUBJECT/surf/%s.$SURF' % pipeline._str_hemi)
                                misc.file_writeOnce('pval_N_statFunc.freeview.txt', '#!/bin/bash\n\nSUBJNAME=$1; SURF=$2\nexport SUBJECT=$SUBJECTS_DIR/$SUBJNAME\nfreeview --surface $SUBJECT/surf/%s.$SURF' % pipeline._str_hemi)
                                misc.file_writeOnce('all.freeview.txt', '#!/bin/bash\n\nSUBJNAME=$1; SURF=$2\nexport SUBJECT=$SUBJECTS_DIR/$SUBJNAME\nfreeview --surface $SUBJECT/surf/%s.$SURF' % pipeline._str_hemi)
                                for pvalROI in l_pval:
                                    if pvalROI != 'entire':
                                        intensity = 127+labelCount*colorMult
                                        if intensity > 255: intensity = 255
                                        misc.file_writeOnce('pval-only.freeview.txt',
                                        ':label=$SUBJECT/label/%s.%s.label:label_color=%d,0,0' % (pipeline._str_hemi, pvalROI, intensity),
                                        mode='a')
                                        misc.file_writeOnce('pval-only.tcl',
                                        'labl_load %s.%s.label ; labl_set_color %d %d 0 0\n' % (pipeline._str_hemi, pvalROI, labelCount, intensity),
                                        mode='a')
                                        misc.file_writeOnce('all.freeview.txt',
                                        ':label=$SUBJECT/label/%s.%s.label:label_color=%d,0,0' % (pipeline._str_hemi, pvalROI, intensity),
                                        mode='a')
                                        misc.file_writeOnce('all.tcl',
                                        'labl_load %s.%s.label ; labl_set_color %d %d 0 0\n' % (pipeline._str_hemi, pvalROI, labelCount, intensity),
                                        mode='a')
                                        labelCount      += 1
                                rLabelCount  = labelCount
                                misc.file_writeOnce('pval-only.tcl',
                                                    tcl_append(pipeline.namespec(), '%s-pval-only' % str_curvFunc),
                                                    mode='a')
                                for statFuncROI in l_statFunc:
                                    if statFuncROI != 'entire':
                                        intensity = 127+(labelCount-rLabelCount)*colorMult
                                        if intensity > 255: intensity = 255
                                        misc.file_writeOnce('statFunc-only.freeview.txt',
                                        ':label=$SUBJECT/label/%s.%s.label:label_color=0,0,%d' % (pipeline._str_hemi, statFuncROI, intensity),
                                        mode='a')
                                        misc.file_writeOnce('statFunc-only.tcl',
                                        'labl_load %s.%s.label ; labl_set_color %d 0 0 %d\n' % (pipeline._str_hemi, statFuncROI, labelCount, intensity),
                                        mode='a')
                                        misc.file_writeOnce('all.freeview.txt',
                                        ':label=$SUBJECT/label/%s.%s.label:label_color=0,0,%d' % (pipeline._str_hemi, statFuncROI, intensity),
                                        mode='a')
                                        misc.file_writeOnce('all.tcl',
                                        'labl_load %s.%s.label ; labl_set_color %d 0 0 %d\n' % (pipeline._str_hemi, statFuncROI, labelCount, intensity),
                                        mode='a')
                                        labelCount += 1
                                rLabelCount  = labelCount
                                misc.file_writeOnce('statFunc-only.tcl',
                                                    tcl_append(pipeline.namespec(), '%s-statFunc-only' % str_curvFunc),
                                                    mode='a')
                                for pNsROI in l_pNs:
                                    if pNsROI != 'entire':
                                        intensity = 127+(labelCount-rLabelCount)*colorMult
                                        if intensity > 255: intensity = 255
                                        misc.file_writeOnce('pval_N_statFunc.freeview.txt',
                                        ':label=$SUBJECT/label/%s.%s.label:label_color=0,%d,0' % (pipeline._str_hemi, pNsROI, intensity),
                                        mode='a')
                                        misc.file_writeOnce('pval_N_statFunc.tcl',
                                        'labl_load %s.%s.label ; labl_set_color %d 0 %d 0\n' % (pipeline._str_hemi, pNsROI, labelCount, intensity),
                                        mode='a')
                                        misc.file_writeOnce('all.freeview.txt',
                                        ':label=$SUBJECT/label/%s.%s.label:label_color=0,%d,0' % (pipeline._str_hemi, pNsROI, intensity),
                                        mode='a')
                                        misc.file_writeOnce('all.tcl',
                                        'labl_load %s.%s.label ; labl_set_color %d 0 %d 0\n' % (pipeline._str_hemi, pNsROI, labelCount, intensity),
                                        mode='a')
                                        labelCount += 1
                                misc.file_writeOnce('pval_N_statFunc.tcl',
                                                    tcl_append(pipeline.namespec(), '%s-pval_N_statFunc' % str_curvFunc),
                                                    mode='a')
                                misc.file_writeOnce('all.tcl',
                                                    tcl_append(pipeline.namespec(), '%s-all' % str_curvFunc),
                                                    mode='a')

                                # Now run the actual tcl files to generate the tiffs.
                                log("Current dir = %s\n" % os.getcwd())
                                str_scriptDir = '/neuro/users/rudolphpienaar/projects/dyslexia-curv-analysis-2/sh'
                                str_execCmd = 'cd %s;  %s/%s -D %s -h %s' % \
                                        (
                                            os.getcwd(),
                                            str_scriptDir,
#                                            pipeline.workingDir(),
                                            "./tksurfer-run.bash",
                                            os.getcwd(),
                                            pipeline._str_hemi,
                                        )
                                log("Shell command = %s\n" % str_execCmd)
                                shellCluster.waitForChild(True)
                                shellCluster.echoStdOut(True)
                                shellCluster.echoStdErr(True)
                                shellCluster(str_execCmd)
        shellCluster.blockOnChild()


        stage.exitCode(0)
        return True
    stage4.def_stage(f_stage4callback, roi=args.l_roi, obj=stage4, pipe=roitag)
    # stage4.def_postconditions(f_blockOnScheduledJobs, obj=stage2,
    #                           blockProcess    = 'mris_label_calc')


    roitaglog = roitag.log()
    roitaglog('INIT: %s\n' % ' '.join(sys.argv))
    roitag.stage_add(stage0)
    roitag.stage_add(stage1)
    roitag.stage_add(stage2)
    roitag.stage_add(stage3)
    roitag.stage_add(stage4)
    roitag.initialize()

    roitag.run()
