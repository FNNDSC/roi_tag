#!/usr/bin/env python

'''

    'rsample.py' simply "subsamples" an experimental directory to generate
    new sub-populations. These can be tested/analyzed in the same manner
    as the parent population to study sample drift.

'''


import  os
import  sys
import  string
import  argparse
import  random
import  tempfile, shutil, glob
from    _common import systemMisc       as misc
from    _common import crun

import  error
import  message
import  stage
import  fnndsc  as base

class FNNDSC_RSAMPLE(base.FNNDSC):
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
        'outDirNotFound':  {
            'action'        : 'attempting to access the <outDir>, ',
            'error'         : 'directory was not found. It will be created.',
            'exitCode'      : 16},    }


    def threshold(self, *args):
        if len(args):
            self._f_threshold   = args[0]
        else:
            return self._f_threshold

    def topDir(self, *args):
        if len(args):
            self._topDir    = args[0]
        else:
            return self._topDir

    def outDir(self, *args):
        if len(args):
            self._outDir    = args[0]
        else:
            return self._outDir

    def annotationDir(self, *args):
        if len(args):
            self._str_annotationDir = args[0]
        else:
            return self._str_annotationDir

    def dirSpecAnnot(self, **kwargs):
        '''
        Return the path to the annotation root directory.

        The 'kwargs' specifies either 'orig' or 'rsampled'
        as the base directory.

        '''
        str_base    = 'orig'
        for key, value in kwargs.iteritems():
            if key == 'base':       str_base    = value
        str_root    = self.topDir()
        if str_base == 'orig':      str_root    = self.topDir()
        if str_base == 'rsampled':  str_root    = self.outDir()
        return '%s/groupCurvAnalysis/%s' % (
                str_root,
                self._str_annotationDir
            )

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


    def subj(self):
        return self._str_subj

    def l_subj(self, *args):
        if len(args):
            self._l_subj = args[0]
        else:
            return self._l_subj

    def l_subjResampled(self, *args):
        if len(args):
            self._l_subjResampled = args[0]
        else:
            return self._l_subjResampled

    def l_ROI(self, *args):
        if len(args):
            self._l_ROI = args[0]
        else:
            return self._l_ROI

    def l_centroidTXT(self, *args):
        if len(args):
            self._l_centroidTXT = args[0]
        else:
            return self._l_centroidTXT

    def __init__(self, **kwargs):
        '''
        Basic constructor. Checks on named input args, checks that files
        exist and creates directories.

        '''
        base.FNNDSC.__init__(self, **kwargs)

        self._lw                        = 120
        self._rw                        = 20
        self._l_subj                    = []
        self._l_subjResampled           = []
        self._l_ROI                     = []
        self._l_centroidTXT             = []

        self._outDir                    = ''
        self._stageslist                = '01'
        self._str_annotationDir         = ''

        # Internal tracking vars
        self._topDir                    = ''
        self._str_subj                  = ''
        self._f_threshold               = 0.5
        self._str_ROI                   = ''
        self._str_centroidTXT           = ''

        for key, value in kwargs.iteritems():
            if key == 'subjectList':    self._l_subj            = value
            if key == 'outDir':         self._outDir            = value
            if key == 'stages':         self._stageslist        = value
            if key == 'threshold':      self._f_threshold       = float(value)
            if key == 'annotation':     self._str_annotationDir = value


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
        os.chdir(self.topDir())
        self.topDir(os.path.abspath(self.topDir()))
        if os.path.isdir(self.outDir()):
            os.chdir(self.outDir())
            self.outDir(os.path.abspath(self.outDir()))
        else:
            error.warn(self, 'outDirNotFound')
            misc.mkdir(self.outDir())


    def run(self):
        '''
        The main 'engine' of the class.

        '''
        base.FNNDSC.run(self)


def synopsis(ab_shortOnly = False):
    scriptName = os.path.basename(sys.argv[0])
    shortSynopsis =  '''
    SYNOPSIS

            %s                                      \\
                            [--stages <stages>]             \\
                            [-o|--outDir <outputRootDir>]   \\
                            [-a|--annotation <annotation>]  \\
                            [-v|--verbosity <verboseLevel>] \\
                            [-s|--stages <stages>]          \\
                            [-r|--random <threshold>]       \\
                            <dir1> <dir2> ... <dir_N>
    ''' % scriptName

    description =  '''
    DESCRIPTION

        `%s' simply copies (or links) passed directories pending
        a random threshold set at --random <threshold>, i.e. a
        random value between 0.0 and 1.0 is selected and if larger
        than <threshold> the directory will be copied (linked) to
        the new <outputRootDir>.

    ARGS

       --stages <stages>
       The stages to execute. This is specified in a string, such as '1234'
       which would imply stages 1, 2, 3, and 4.

       The special keyword 'all' can be used to turn on all stages.

       --outDir <outputRootDir>
       The output root directory containing the re-sampled experiment.

       --annotation <annotation>
       The original annotation directory to resample.

       --random <threshold>
       The pval cutoffs to consider. In practice, this is always 'le1,le5'

       <dir1> <dir2> ... <dirN>
       The directories (subjects) to process.

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

    parser.add_argument('l_subj',
                        metavar='SUBJs', nargs='+',
                        help='Subject directories to process')
    parser.add_argument('--verbosity', '-v',
                        dest='verbosity',
                        action='store',
                        default=0,
                        help='verbosity level')
    parser.add_argument('--output', '-o',
                        dest='outDir',
                        action='store',
                        default='rsample',
                        help='output root directory')
    parser.add_argument('--stages', '-s',
                        dest='stages',
                        action='store',
                        default='01',
                        help='analysis stages')
    parser.add_argument('--random', '-r',
                        dest='threshold',
                        action='store',
                        default='0.50',
                        help='directory selection threshold')
    parser.add_argument('--annot', '-a',
                        dest='annotation',
                        action='store',
                        default='',
                        help='annotation directory to process')

    args = parser.parse_args()

    # A generic "shell"
    OSshell = crun.crun()
    OSshell.echo(False)
    OSshell.echoStdOut(False)
    OSshell.detach(False)
    OSshell.waitForChild(True)

    rsample = FNNDSC_RSAMPLE(
                        subjectList     = args.l_subj,
                        outDir          = args.outDir,
                        stages          = args.stages,
                        threshold       = args.threshold,
                        annotation      = args.annotation,
                        logTo           = 'rsample.log',
                        syslog          = True,
                        logTee          = True)

    rsample.verbosity(args.verbosity)
    pipeline    = rsample.pipeline()
    pipeline.poststdout(True)
    pipeline.poststderr(True)

    rsample.topDir(os.getcwd())

    stage0 = stage.Stage(
                        name            = 'rsample-0-filter',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'rsample-0-filter.log',
                        logTee          = True,
                        )
    def f_stage0callback(**kwargs):
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        log             = stage.log()
        lst_subj        = pipeline.l_subj()

        if os.path.isdir(pipeline.outDir()):
            log('Existing outDir tree found... deleting...\n')
            shutil.rmtree(pipeline.outDir())
        log('Creating output directory...\n')
        OSshell('mkdir -p %s' % pipeline.outDir())
        print(OSshell.stdout())
        if OSshell.exitCode() != 0: error.fatal(pipeline, 'outDirNotCreate')

        os.chdir(pipeline.outDir())
        pipeline.outDir(os.getcwd())

        passcount   = 0
        failcount   = 0
        total       = 0
        for pipeline._str_subj  in lst_subj:
            f_cutoff    = random.random()
            log('thresholding %s...' % pipeline._str_subj)
            if f_cutoff < pipeline.threshold():
                OSshell('ln -s %s/%s .' % (
                        pipeline.topDir(), pipeline._str_subj))
                log('[ passed ]\n', rw=20, syslog=False)
                passcount += 1
            else:
                log('[ failed ]\n', rw=20, syslog=False)
                failcount += 1
            total += 1
        f_passPerc  = float(passcount)/float(total)*100
        f_failPerc  = float(failcount)/float(total)*100
        log('passed: %d (%5.2f%s)\n' % (passcount, f_passPerc, '%'))
        log('failed: %d (%5.2f%s)\n' % (failcount, f_failPerc, '%'))
        stage.exitCode(0)
        return True
    stage0.def_stage(f_stage0callback, obj=stage0, pipe=rsample)
#    stage0.def_postconditions(f_stageShellExitCode, obj=stage0)

    stage1 = stage.Stage(
                        name            = 'rsample-1-ROIcopy',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'rsample-1-ROIcopy.log',
                        logTee          = True,
                        )
    def f_stage1callback(**kwargs):
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        log             = stage.log()
        lst_subj        = pipeline.l_subj()

        os.chdir(pipeline.outDir())
        misc.mkdir('groupCurvAnalysis')
        os.chdir('groupCurvAnalysis')
        misc.mkdir(pipeline.annotationDir())
        os.chdir(pipeline.dirSpecAnnot(base='orig'))
        for bashFile in glob.glob(r'*.bash'):
            shutil.copy(bashFile, pipeline.dirSpecAnnot(base='rsampled'))
        for lstFile in glob.glob(r'*.lst'):
            shutil.copy(lstFile, pipeline.dirSpecAnnot(base='rsampled'))

        str_CMDoriginalROI  = "find . -maxdepth 1 -mindepth 1 -type d | grep -v ROItag | awk -F\/ '{print $2}'"
        pipeline.l_ROI(OSshell(str_CMDoriginalROI)[0].strip().split('\n'))
        lst_ROI         = pipeline.l_ROI()
        os.chdir(pipeline.dirSpecAnnot(base='rsampled'))
        # print(os.getcwd())

        total = 0
        for pipeline._str_ROI in lst_ROI:
            os.chdir(pipeline.dirSpecAnnot(base='rsampled'))
            log('Creating ROI dir %s\n' % (pipeline._str_ROI))
            misc.mkdir(pipeline._str_ROI)
            os.chdir(pipeline._str_ROI)
            tagCentroidCMD = 'ls %s/%s/*centroid*txt' % (   pipeline.dirSpecAnnot(base='orig'),
                                                            pipeline._str_ROI)
            pipeline.l_centroidTXT(OSshell(tagCentroidCMD)[0].strip().split('\n'))
            for pipeline._str_centroidTXT in pipeline.l_centroidTXT():
                log('\tCopying original sample %s\n' % os.path.basename(pipeline._str_centroidTXT))
                shutil.copy(pipeline._str_centroidTXT,
                            os.path.join(pipeline.outDir(), 'groupCurvAnalysis',
                            pipeline.annotationDir(), pipeline._str_ROI))

            total += 1
        stage.exitCode(0)
        return True
    stage1.def_stage(f_stage1callback, obj=stage1, pipe=rsample)

    stage2 = stage.Stage(
                        name            = 'rsample-2-centroidFilter',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = 'rsample-2-centroidFilter.log',
                        logTee          = True,
                        )
    def f_stage2callback(**kwargs):
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val
        log             = stage.log()
        lst_subj        = pipeline.l_subj()

        os.chdir(pipeline.outDir())
        str_CMDsubj     = "find . -maxdepth 1 -mindepth 1 -type l | awk -F\/ '{print $2}' | sort "
        pipeline.l_subjResampled(OSshell(str_CMDsubj)[0].strip().split('\n'))
        lst_subjResampled   = pipeline.l_subjResampled()
        lst_subjResampled.insert(0, 'Subj')
        os.chdir(pipeline.dirSpecAnnot(base='rsampled'))
        str_CMDoriginalROI  = "find . -maxdepth 1 -mindepth 1 -type d | grep -v ROItag | awk -F\/ '{print $2}'"
        pipeline.l_ROI(OSshell(str_CMDoriginalROI)[0].strip().split('\n'))
        lst_ROI         = pipeline.l_ROI()

        total = 0
        for pipeline._str_ROI in lst_ROI:
            os.chdir('%s/%s' % (pipeline.dirSpecAnnot(base='rsampled'), pipeline._str_ROI))
            str_CMDcentroid = "ls -1 *centroid*txt"
            pipeline.l_centroidTXT(OSshell(str_CMDcentroid)[0].strip().split('\n'))
            log('Switching to ROI %s\n' % pipeline._str_ROI)
            for pipeline._str_centroidTXT in pipeline.l_centroidTXT():
                log('\tProcessing %s\n' % pipeline._str_centroidTXT)
                l_meas  = [line.strip() for line in open(pipeline._str_centroidTXT, 'r')]
                filt    = [x for x in l_meas if x.split()[0] in lst_subjResampled]
                FILEfilt= open('%s.filt' % pipeline._str_centroidTXT, 'w')
                for line in filt:
                    FILEfilt.write('        %s\n' % line)
                # log('Saved in %s\n' % os.getcwd())
                FILEfilt.close()
                shutil.move('%s.filt' % pipeline._str_centroidTXT, pipeline._str_centroidTXT)
            total += 1
        stage.exitCode(0)
        return True
    stage2.def_stage(f_stage2callback, obj=stage2, pipe=rsample)

    rsamplelog = rsample.log()
    rsamplelog('INIT: %s\n' % ' '.join(sys.argv))
    rsample.stage_add(stage0)
    rsample.stage_add(stage1)
    rsample.stage_add(stage2)
    rsample.initialize()

    rsample.run()
