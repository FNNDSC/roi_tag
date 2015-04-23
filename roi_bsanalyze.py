#!/usr/bin/env python

'''

    'roi_bsanalyze.py' analyzes the output directories of a bootstrapped
    curvature analysis run subsample set and generates an ROItag tree with
    each leaf node containing the probability of ROI occurence based on
    observations across all the bootstraps.

'''


import  os
import  sys
import  string
import  argparse
import  tempfile, shutil
import  json
from    json    import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.2f')

from    _common import systemMisc       as misc
from    _common import crun

import  error
import  message
import  stage
import  fnndsc  as base
import  C_snode


class FNNDSC_ROIBSANALYZE(base.FNNDSC):
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
        str_curvFunc         = ','.join(self.l_curvFunc())
        str_dirSpecSubsample = '%s/groupCurvAnalysis/%s/ROItag-%s/%s/%s/%s/%s/%s/%s' % (
                                        self._d_subsampleDir[self._str_subsample],
                                        self._str_annotation,
                                        str_curvFunc,
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
        self._stageslist                = '12'

        # Internal tracking vars
        self._str_pval                  = ''
        self._str_group                 = ''
        self._str_roi                   = ''
        self._str_hemi                  = ''
        self._str_surface               = ''
        self._str_statFunc              = ''
        self._str_subsample             = ''
        self._str_annotation            = ''
        self._str_ctype                 = ''

        self._topDir                    = ''

        self._dtree                     = C_snode.C_stree(['/'])

        for key, value in kwargs.iteritems():
            if key == 'rsampleList':    self._l_subsample       = value
            if key == 'annot':          self._l_annotation      = value.split(',')
            if key == 'outDir':         self._outDir            = value
            if key == 'stages':         self._stageslist        = value
            if key == 'pval':           self._l_pval            = value.split(',')
            if key == 'group':          self._l_group           = value.split(',')
            if key == 'surface':        self._l_surface         = value.split(',')
            if key == 'statFunc':       self._l_statFunc        = value.split(',')
            if key == 'curvFunc':       self._l_curvFunc        = value.split(',')
            if key == 'hemi':           self._l_hemi            = value.split(',')


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

        # Set absolute dir specs for topDir and outDir
        log = self.log()
        log('Setting dir specs...\n')
        os.chdir(self.topDir())
        self.topDir(os.path.abspath(self.topDir()))
        if os.path.isdir(self.outDir()):
            log('Existing output directory found. Deleting and recreating...\n')
            shutil.rmtree(self.outDir())
        misc.mkdir(self.outDir())
        os.chdir(self.outDir())
        self.outDir(os.getcwd())
        os.chdir(self.topDir())

        # check dirs for each subsample dir
        # print(os.getcwd())
        # print(self._l_subsample)
        for directory in self._l_subsample:
            if os.path.isdir(os.path.join(self.topDir(), directory)):
                os.chdir(directory)
                self._l_subsampleDir.append(os.path.abspath(os.getcwd()))
                os.chdir(self.topDir())
            else:
                error.fatal(self, 'subsampleDirNotFound')
        self._d_subsampleDir = dict(zip(self._l_subsample, self._l_subsampleDir))
        # print(self._d_subsampleDir)

        self.dtree_build()
        self._dtree.tree_metaData_print(False)
        # print(self._dtree)

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
                            [-o|--outDir <outputRootDir>]   \\
                            [-a|--annot <annotation>]       \\
                            [-v|--verbosity <verboseLevel>] \\
                            [-s|--stages <stages>]          \\
                            [-p|--pval <pvalCutoffList>]    \\
                            [-g|--group <groupList>]        \\
                            [-S|--surface <surfaceList>]    \\
                            [-f|--statFunc <statFuncList>]  \\
                            [-c|--curvFunc <curvFuncList>]  \\
                            [-m|--hemi <hemiList>]          \\
                            <subsample1> <subsample2> ... <subsampleN>
    ''' % scriptName

    description =  '''
    DESCRIPTION

        `%s' analyzes trees of subsampled subject spaces in a bootstrap type
        approach.

        For each subsample tree, the ROI tags at each leaf node for a single
        branch descent are counted, and the probability of any single ROI tag
        across all subsamples is calculated.

        This tree structure is identical to the ROItag tree:

        <group>
          |
          +--<pval>
            |
            +--<statFunc>
              |
              +--<surface>
                |
                +--<hemi>
                  |
                  +-- curv
                  | +-- pval
                  | +-- statFunc
                  | +-- pval-statFunc
                  +-- thickness
                    +-- pval
                    +-- statFunc
                    +-- pval-statFunc

        Each of the above subdirectories will contain a list of label ROI files
        that satisfy the named constraint.

    ARGS

       --outDir <outputRootDir>
       The root directory to contain the results of the bootstrap.

       --stages <stages>
       The stages to execute. This is specified in a string, such as '1234'
       which would imply stages 1, 2, 3, and 4.

       The special keyword 'all' can be used to turn on all stages.

       --annot <annotation>
       The annotation to analyze. Assumes that this directory contains a subtree,
       'ROItag-<curvFunc>'.

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

       <subsample1> <subsample2> ... <subsampleN>
       The subsamples to process.

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

    parser.add_argument('l_resample',
                        metavar='Samples', nargs='+',
                        help='resampled experiment list to process')
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
                        default='ROIbsanalyze',
                        help='output root directory')
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

    args = parser.parse_args()

    # A generic "shell"
    OSshell = crun.crun()
    OSshell.echo(False)
    OSshell.echoStdOut(False)
    OSshell.detach(False)
    OSshell.waitForChild(True)

    roibsanalyze = FNNDSC_ROIBSANALYZE(
                        rsampleList     = args.l_resample,
                        annot           = args.annotation,
                        outDir          = args.outDir,
                        stages          = args.stages,
                        pval            = args.pval,
                        group           = args.group,
                        surface         = args.surface,
                        statFunc        = args.statFunc,
                        curvFunc        = args.curvFunc,
                        hemi            = args.hemi,
                        logTo           = 'ROIbsanalyze.log',
                        syslog          = True,
                        logTee          = True)

    roibsanalyze.verbosity(args.verbosity)
    pipeline    = roibsanalyze.pipeline()
    pipeline.poststdout(True)
    pipeline.poststderr(True)
    roibsanalyze.topDir(os.getcwd())

    stage0 = stage.Stage(
                        name            = 'ROIbsanalyze-0-init',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIbsanalyze-0-init.log',
                        logTee          = True,
                        )
    def f_stage0callback(**kwargs):
        str_cwd         =  os.getcwd()
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        log             = stage.log()
        lst_annotation  = pipeline.l_annotation()
        lst_pval        = pipeline.l_pval()
        lst_group       = pipeline.l_group()
        lst_hemi        = pipeline.l_hemisphere()
        lst_surface     = pipeline.l_surface()
        lst_statFunc    = pipeline.l_statFunc()
        lst_ctype       = ['thickness', 'curv']
        lst_leaf        = ['pval', 'statFunc', 'pval-only', 'statFunc-only', 'pval_N_statFunc']

        for pipeline._str_annotation in lst_annotation:
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

    stage0.def_stage(f_stage0callback, obj=stage0, pipe=roibsanalyze)

    stage1 = stage.Stage(
                        name            = 'ROIbsanalyze-1-probability',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'ROIbsanalyze-1-probability.log',
                        logTee          = True
                        )
    stage1.def_preconditions(lambda **x: True)
    def f_stage1callback(**kwargs):
        ''' Sorts intersect/xor of classes '''
        str_cwd         =  os.getcwd()
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        log             = stage.log()
        lst_annotation  = pipeline.l_annotation()
        lst_subsample   = pipeline.l_subsample()
        lst_pval        = pipeline.l_pval()
        lst_group       = pipeline.l_group()
        lst_hemi        = pipeline.l_hemisphere()
        lst_surface     = pipeline.l_surface()
        lst_statFunc    = pipeline.l_statFunc()
        lst_ctype       = pipeline.l_ctype()

        c               = pipeline.dtree()

        for pipeline._str_annotation in lst_annotation:
            for pipeline._str_group in lst_group:
                for pipeline._str_pval in lst_pval:
                    for pipeline._str_statFunc in lst_statFunc:
                        for pipeline._str_surface in lst_surface:
                            for pipeline._str_hemi in lst_hemi:
                                for pipeline._str_ctype in lst_ctype:
                                    l_sum       = []
                                    d_occurence = {}
                                    for pipeline._str_subsample in lst_subsample:
                                        log('Processing %s %s for %s...' % (pipeline.dirSpecPartial(),
                                                                        pipeline._str_ctype,
                                                                        pipeline._str_subsample), lw=170)
                                        #print(pipeline.dirSpecSubsample())
                                        os.chdir(pipeline.dirSpecSubsample())
                                        l_roi       = OSshell("cat all.tcl | grep label | awk '{print $2}'")[0].strip().split('\n')
                                        log('[ %02d ]\n' % len(l_roi), syslog=False, rw=20)
                                        #print(l_roi)
                                        c.cdnode(pipeline.dtreeDir())
                                        c.touch('l_roi', l_roi)
                                        c.touch('path', pipeline.dtreeDir())
                                        l_sum += l_roi
                                    c.cdnode('../')
                                    c.touch('l_sum', l_sum)
                                    l_uniq = sorted(set(l_sum))
                                    for roi in l_uniq:
                                        d_occurence[roi] = float(l_sum.count(roi))/float(len(lst_subsample))*100
                                    c.touch('d_occurence', d_occurence)
                                    os.chdir(pipeline.d2path())
                                    # print(os.getcwd())
                                    # print(json.dumps(d_occurence, sort_keys=True, indent=4,
                                    #                 separators=(',',': ')))
                                    json.dump(d_occurence, open('%s/occurence.txt' % os.getcwd(), "w"),
                                                sort_keys=True,
                                                indent=4,
                                                separators=(',',': ')
                                                )
        stage.exitCode(0)
        return True
    stage1.def_stage(f_stage1callback, obj=stage1, pipe=roibsanalyze)


    roibsanalyzelog = roibsanalyze.log()
    roibsanalyzelog('INIT: %s\n' % ' '.join(sys.argv))
    roibsanalyze.stage_add(stage0)
    roibsanalyze.stage_add(stage1)
    roibsanalyze.initialize()

    roibsanalyze.run()
