#!/usr/bin/env python

'''

    'roi_gcibs.py' compares two groups informed by an a priori bootstrap analysis.

'''


import  os
import  sys
import  argparse
import  tempfile, shutil
import  json
import  pprint
import  copy

from    collections     import  defaultdict


from    _common         import  systemMisc      as misc
from    _common         import  crun

import  error
import  message
import  stage
import  fnndsc  as base

class FNNDSC_roigcibs(base.FNNDSC):
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
        """
        Return the dirSpec based on internal pipeline._str_* variables
        """
        return '%s/%s/%s/%s/%s/%s/%s' % (
                                        self.outDir(),
                                        self._str_annotation,
                                        self._str_group,
                                        self._str_pval,
                                        self._str_statFunc,
                                        self._str_surface,
                                        self._str_hemi
        )

    def dirSpecPartial(self):
        """
        Return the dirSpec based on internal pipeline._str_* variables w/o
        the leading directories.
        """
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
        return '%s%s%s%s%s%s%s%s%s%s%s' % (
                                        self._str_annotation,   str_sep,
                                        self._str_group,        str_sep,
                                        self._str_pval,         str_sep,
                                        self._str_statFunc,     str_sep,
                                        self._str_surface,      str_sep,
                                        self._str_hemi
        )

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
        """
        Basic constructor. Checks on named input args, checks that files
        exist and creates directories.

        """
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
        self._l_annot                   = []

        self._outDir                    = ''
        self._workingDir                = ''
        self._stageslist                = '12'

        self._f_lowerBoundHard          = 0.0
        self._f_lowerBoundSoft          = 0.0
        self._f_upperBoundSoft          = 0.0

        # Internal tracking vars
        self._str_pval                  = ''
        self._str_group                 = ''
        self._str_roi                   = ''
        self._str_hemi                  = ''
        self._str_surface               = ''
        self._str_statFunc              = ''
        self._str_curvFunc              = ''
        self._str_annotation            = ''

        self._topDir                    = ''

        self._d_bootstrapOccurrence     = Tree()
        self._d_bootstrapThreshold      = Tree()
        self._d_bootstrapFiltered       = Tree()

        # Scheduler std out/err dirs
        self._str_schedulerStdOutDir    = '~/scratch'
        self._str_schedulerStdErrDir    = '~/scratch'

        self._b_clobber                 = False

        for key, value in kwargs.iteritems():
            if key == 'outDir':                 self._outDir                    = value
            if key == 'workingDir':             self._workingDir                = value
            if key == 'stages':                 self._stageslist                = value
            if key == 'curvFunc':               self._l_curvFunc                = value.split(':')
            if key == 'pval':                   self._l_pval                    = value.split(',')
            if key == 'group':                  self._l_group                   = value.split(',')
            if key == 'surface':                self._l_surface                 = value.split(',')
            if key == 'statFunc':               self._l_statFunc                = value.split(',')
            if key == 'hemi':                   self._l_hemi                    = value.split(',')
            if key == 'annot':                  self._l_annot                   = value.split(',')
            if key == 'lowerBoundSoft':         self._f_lowerBoundSoft          = float(value)
            if key == 'lowerBoundHard':         self._f_lowerBoundHard          = float(value)
            if key == 'upperBoundSoft':         self._f_upperBoundSoft          = float(value)
            if key == 'schedulerStdOutDir':     self._str_schedulerStdOutDir    = value
            if key == 'schedulerStdErrDir':     self._str_schedulerStdErrDir    = value

        if not os.path.isdir(self._workingDir): errorFatal(self, 'workingDirNotExist')


    def initialize(self):
        """
        This method provides some "post-constructor" initialization. It is
        typically called after the constructor and after other class flags
        have been set (or reset).

        """

        # Set the stages
        self._pipeline.stages_canRun(False)
        lst_stages = list(self._stageslist)
        for index in lst_stages:
            stage = self._pipeline.stage_get(int(index))
            stage.canRun(True)

    def run(self):
        """
        The main 'engine' of the class.

        """
        base.FNNDSC.run(self)

    def innerLoop(self, func_callBack, *args, **callBackArgs):
        '''

        A loop function that calls func_callBack(**callBackArgs)
        at the innermost loop the nested data dictionary structure.

        The loop order:

            annotation, group, pval, statFunc, surface, hemi

        Note that internal tracking object variables, _str_gid ... _str_ctype
        are automatically updated by this method.

        The **callBackArgs is a generic dictionary holder that is interpreted
        by both this loop controller and also passed down to the callback
        function.

        In the context of the loop controller, loop conditions can
        be changed by passing appropriately name args in the
        **callBackArgs structure.

        '''
        ret             = True
        _str_log        = ''

        for key, val in callBackArgs.iteritems():
            if key == 'hemi':           self._l_hemi        = val
            if key == 'surface':        self._l_surface     = val
            if key == 'curv':           self._l_curvFunc    = val
            if key == 'group':          self._l_group       = val
            if key == 'log':            _str_log            = val

        if len(_str_log): self._log(_str_log)

        for self._str_annotation in self._l_annot:
            for self._str_group in self._l_group:
                for self._str_pval in self._l_pval:
                    for self._str_statFunc in self._l_statFunc:
                        for self._str_surface in self._l_surface:
                            for self._str_hemi in self._l_hemi:
                                for self._str_curvFunc in self._l_curvFunc:
                                    ret = func_callBack(**callBackArgs)

        if len(_str_log): self._log('[ ok ]\n', rw=self._rw)
        return ret

    def outputDirTree_build(self, **kwargs):
        '''Build the tree structure containing output images
        '''
        OSshell('mkdir -p %s/%s' % (
                self.dirSpec(),
                self._str_curvFunc
        ))

    def tcl_append(self, str_prefix, str_suffix, str_hemi):
        """
        Append some text to each tcl file to save snapshots of different brain aspects.
        """
        if str_hemi == "lh":
            frontal_or_distal = "frontal"
            distal_or_frontal = "distal"
            medial_or_lateral = "medial"
            lateral_or_medial = "lateral"
        if str_hemi == "rh":
            frontal_or_distal = "distal"
            distal_or_frontal = "frontal"
            medial_or_lateral = "lateral"
            lateral_or_medial = "medial"
        str_content = '''
        # Initial lateral view
        read_binary_curv
        redraw
        save_tiff %s/lateral-%s.tiff

        # inferior view
        rotate_brain_x 90
        redraw
        save_tiff %s/inferior-%s.tiff

        # superior view
        rotate_brain_x -180
        redraw
        save_tiff %s/superior-%s.tiff

        # reset
        rotate_brain_x 90

        # %s view
        rotate_brain_y 90
        redraw
        save_tiff %s/%s-%s.tiff

        # medial view
        rotate_brain_y 90
        redraw
        save_tiff %s/medial-%s.tiff

        # %s view
        rotate_brain_y 90
        redraw
        save_tiff %s/%s-%s.tiff

        exit 0
        ''' % ( str_prefix, str_suffix,
                str_prefix, str_suffix,
                str_prefix, str_suffix,
                distal_or_frontal,
                str_prefix, distal_or_frontal, str_suffix,
                str_prefix, str_suffix,
                frontal_or_distal,
                str_prefix, frontal_or_distal, str_suffix)
        return str_content

    def static_vars(**kwargs):
        def decorate(func):
            for k in kwargs:
                setattr(func, k, kwargs[k])
            return func
        return decorate

    def labelScript_process(self, **kwargs):
        """
        Write the tcl files to display filtered ROIs
        :param kwargs:
        :return:
        """
        spec    = self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc].keys()[0]
        innerDict = self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec]
        # print(innerDict)
        str_dirSpec         = '%s/%s'       % (self.dirSpec(), self._str_curvFunc)
        os.chdir(str_dirSpec)
        str_fileStem        = "%s-%s"       % (self.namespec("-"), self._str_curvFunc)
        str_TCLfileName     = '%s.tcl'      % (str_fileStem)
        str_JSONfileName    = '%s.json'     % (str_fileStem)
        self._log("\n")
        str_currentDir      = os.getcwd()
        l_currentDir        = str_currentDir.split('/')
        l_workingDir        = l_currentDir[-8:-1]
        l_workingDir.append(l_currentDir[-1])
        index               = 0
        self._log("Current dir:    %s\n" % '/'.join(l_workingDir))
        self._log('Creating tcl file:   %s...\n' % str_TCLfileName)
        with open(str_JSONfileName, 'w') as JSONfile:
            json.dump(innerDict, JSONfile, indent=4, sort_keys=True)
            self._log('Creating JSON file:  %s...\n' % str_JSONfileName)
        for key,val in innerDict.iteritems():
            if val <= self._f_lowerBoundSoft:
                misc.file_writeOnce(str_TCLfileName,
                                    'labl_load %s ; labl_set_color %d 0 0 %d\n' %\
                                    (key, index, 2*int(val)), mode='a')
            if val > self._f_lowerBoundSoft and val < self._f_upperBoundSoft:
                misc.file_writeOnce(str_TCLfileName,
                                    'labl_load %s ; labl_set_color %d 0 %d 0\n' %\
                                    (key, index, 2*int(val), 2*int(val)), mode='a')
            if val >= self._f_upperBoundSoft:
                misc.file_writeOnce(str_TCLfileName,
                                    'labl_load %s ; labl_set_color %d %d 0 0\n' %\
                                    (key, index, 2*int(val)), mode='a')
            index += 1
        misc.file_writeOnce(str_TCLfileName, self.tcl_append(str_dirSpec, str_fileStem, self._str_hemi), mode='a')
        str_scriptDir   = '/neuro/users/rudolphpienaar/projects/dyslexia-curv-analysis-2/sh'
        str_subjectDir  = '/neuro/users/rudolphpienaar/projects/dyslexia-curv-analysis-2/results/6-exp-dyslexia-run'
        str_execCmd = 'cd %s;  %s/%s -S %s -D %s -h %s' % \
                      (
                          os.getcwd(),
                          str_scriptDir,
                          "./tksurfer-run.bash",
                          str_subjectDir,
                          os.getcwd(),
                          self._str_hemi,
                      )
        self._log("Shell command = %s\n" % str_execCmd)
        self._log('Executing tcl file...\n')
        OSshell(str_execCmd)
        return {"return": "ok"}

    def bootstrap_filteredDictionaryBuild(self, **kwargs):
        '''Filters group compared results.

        :param kwargs:
        :return:
        '''
        spec    = self._d_bootstrapOccurrence[self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc].keys()[0]
        self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc] = self._d_bootstrapOccurrence[self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc].copy()
        for key in self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec].keys():
            if self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] < \
                self._d_bootstrapThreshold['%s-threshold' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] or \
                self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] <= self._f_lowerBoundHard:
                self._d_bootstrapFiltered['%s-filtered' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec].pop(key, None)
        return {"return": self._d_bootstrapFiltered}

    def bootstrap_thresholdDictionaryBuild(self, **kwargs):
        '''Sum the intra-group occurrences for a lower confidence bound

        :param kwargs:
        :return:
        '''
        str_thresholdOperation      = "sum"
        for kwarg,val in kwargs.iteritems():
            if kwarg == 'threshold':    str_thresholdOperation = val
        str_g1  = self._str_group[0]
        str_g2  = self._str_group[1]
        spec    = self._d_bootstrapOccurrence[str_g1][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc].keys()[0]
        self._d_bootstrapThreshold['%s-threshold' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc] = copy.deepcopy(self._d_bootstrapOccurrence[self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc])
        for key in self._d_bootstrapOccurrence[str_g1][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec].keys():
            if str_thresholdOperation == "sum":
                self._d_bootstrapThreshold['%s-threshold' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] = \
                    self._d_bootstrapOccurrence[str_g1][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] + \
                    self._d_bootstrapOccurrence[str_g2][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key]
            if str_thresholdOperation == "max":
                if self._d_bootstrapOccurrence[str_g1][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] >= self._d_bootstrapOccurrence[str_g2][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key]:
                    self._d_bootstrapThreshold['%s-threshold' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] = \
                        self._d_bootstrapOccurrence[str_g1][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key]
                else:
                    self._d_bootstrapThreshold['%s-threshold' % self._str_group][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key] = \
                        self._d_bootstrapOccurrence[str_g2][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec][key]
        return {"return": self._d_bootstrapThreshold}

    def bootstrap_occurrenceDictionariesBuild(self, **kwargs):
        """Build the occurrence dictionaries:

        This method captures the bootstrap occurrence dictionaries for the
        comparison group as well as the two intra-group variance occurrences.

        """
        str_g1 = self._str_group[0]
        str_g2 = self._str_group[1]
        str_g3 = self._str_group
        l_subgroup = [str_g1, str_g2, str_g3]
        for str_subgroup in l_subgroup:
            os.chdir(self._workingDir)
            if str_subgroup in [str_g1, str_g2]:
                str_g           = "%s000" % str_subgroup
                str_compGroup   = '12'
            else:
                str_g           = '6'
                str_compGroup   = self._str_group
            str_bsDir   = "bootstrap-%s-%s/%s/%s" % (str_g,
                                                    self._str_curvFunc,
                                                    self._str_annotation,
                                                    str_compGroup)
            os.chdir(str_bsDir)
            if self._str_curvFunc == "thickness":
                str_cfilt = "thickness"
            else:
                str_cfilt = "curv"
            self._log('Parsing occurrence data for %2d.%s.%s.%s.%s\n' % (int(str_subgroup), self._str_annotation, self._str_hemi, self._str_surface, self._str_curvFunc))
            OSshell('find . -iname occurence.txt | grep %s | grep %s | grep %s | ../../../../sh/ocfilt.sh -t 0 | python -m json.tool ' % \
                    (str_cfilt, self._str_hemi, self._str_surface))
            self._d_bootstrapOccurrence[str_subgroup][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc] = json.loads(OSshell.stdout())
            spec = self._d_bootstrapOccurrence[str_subgroup][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc].keys()[0]
            self._d_bootstrapOccurrence[str_subgroup][self._str_annotation][self._str_hemi][self._str_surface][self._str_curvFunc][spec].pop("", None)
        return {"return": self._d_bootstrapOccurrence}


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
                            [-a|--annot <annotList>]        \\
                            [-m|--hemi <hemiList>]          \\
                            [--schedulerStdOutDir <dir>]    \\
                            [--schedulerStdErrDir <dir>]    \\

    ''' % scriptName

    description =  '''
    DESCRIPTION

        `%s' performs a group comparison informed by an a priori bootstrap
        analysis. The bootstrap analysis provides a bound on the feature
        variance within a group, and the bottom bound of the threshold
        of significance becomes the sum of the underlying group variances.


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
       The curvature bootstrap analysis to use. Typically
       'K,BE,S', 'H,K1,K2', 'thickness'.

       --hemi <hemiList>
       The hemispheres to process. In practice, this is always 'lh,rh'.

       --annot <annotList>
       The annotation list to process.

       --threshold <thresholdOperation>
       The operation to apply in thresholding variances from the underlying
       bootstrap variances. Either "sum" or "max".

       --lowerBound <lowerBound>
       The lower bound for filtered (and thresholded) comparisons.

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

       For a distributed experiment, delete *all* existing roigcibs trees
       *before* running the experiment!

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

    parser.add_argument('--verbosity', '-v',
                        dest='verbosity',
                        action='store',
                        default=0,
                        help='verbosity level')
    parser.add_argument('--output', '-o',
                        dest='outDir',
                        action='store',
                        default='roigcibs',
                        help='output root directory')
    parser.add_argument('--workingDir', '-w',
                        dest='workingDir',
                        action='store',
                        default='./',
                        help='output working directory')
    parser.add_argument('--clobber', '-C',
                        dest='clobber',
                        action='store_true',
                        default=False,
                        help='if specified, do not erase existing output dir if found.')
    parser.add_argument('--stages', '-s',
                        dest='stages',
                        action='store',
                        default='01',
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
    parser.add_argument('--annot', '-a',
                        dest='annot',
                        action='store',
                        default='aparc.annot',
                        help='comma separated annotation list to process')
    parser.add_argument('--threshold', '-t',
                        dest='threshold',
                        action='store',
                        default='max',
                        help='the threshold operation -- "max" or "sum"')
    parser.add_argument('--lowerBoundHard', '-L',
                        dest='lowerBoundHard',
                        action='store',
                        default=0.0,
                        help='the hard lower bound for filtered occurrences.')
    parser.add_argument('--lowerBoundSoft', '-l',
                        dest='lowerBoundSoft',
                        action='store',
                        default=80.0,
                        help='the soft lower bound for filtered occurrences.')
    parser.add_argument('--upperBoundSoft', '-u',
                        dest='upperBoundSoft',
                        action='store',
                        default=94.0,
                        help='the soft upper bound for filtered occurrences.')
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

    Tree     = lambda:  defaultdict(Tree)

    roigcibs = FNNDSC_roigcibs(
                        outDir                  = '%s' % (args.outDir),
                        workingDir              = args.workingDir,
                        stages                  = args.stages,
                        pval                    = args.pval,
                        group                   = args.group,
                        surface                 = args.surface,
                        curvFunc                = args.curvFunc,
                        statFunc                = args.statFunc,
                        hemi                    = args.hemi,
                        annot                   = args.annot,
                        lowerBoundHard          = args.lowerBoundHard,
                        lowerBoundSoft          = args.lowerBoundSoft,
                        upperBoundSoft          = args.upperBoundSoft,
                        schedulerStdOutDir      = args.schedulerStdOutDir,
                        schedulerStdErrDir      = args.schedulerStdErrDir,
                        logTo                   = '%s/roigcibs.log' % args.workingDir,
                        syslog                  = True,
                        logTee                  = True)

    roigcibs.clobber(args.clobber)
    roigcibs.verbosity(args.verbosity)
    pipeline    = roigcibs.pipeline()
    pipeline.poststdout(True)
    pipeline.poststderr(True)
    os.chdir(roigcibs._workingDir)
    roigcibs._workingDir = os.getcwd()
    roigcibs.topDir(os.getcwd())


    stage0 = stage.Stage(
        name            = 'roigcibs-0-init',
        fatalConditions = True,
        syslog          = True,
        logTo           = '%s/roigcibs-0-init.log' % args.workingDir,
        logTee          = True,
    )
    def f_stage0callback(**kwargs):
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val

        log             = stage._log

        os.chdir(pipeline._workingDir)
        if os.path.isdir(pipeline.outDir()) and not pipeline.clobber():
            log('Existing outDir tree found... deleting...\n')
            shutil.rmtree(pipeline.outDir())
        OSshell('mkdir -p %s' % pipeline.outDir())
        os.chdir(pipeline.outDir())
        pipeline.outDir(os.getcwd())
        if OSshell.exitCode() != 0: error.fatal(pipeline, 'outDirNotCreate')

        d_ret = pipeline.innerLoop(
            pipeline.outputDirTree_build,
            log = "Building output directory tree...\n"
        )
        stage.exitCode(0)
        return True
    stage0.def_stage(f_stage0callback, obj=stage0, pipe=roigcibs)

    stage1 = stage.Stage(
                        name            = 'roigcibs-1-filter',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/roigcibs-1-filter.log' % args.workingDir,
                        logTee          = True,
                        )
    def f_stage1callback(**kwargs):
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val

        os.chdir(pipeline._workingDir)

        d_ret = pipeline.innerLoop(
                    pipeline.bootstrap_occurrenceDictionariesBuild,
                    log = "Parsing bootstrap occurrences...\n"
                )
        d_bootstrapOccurrence = d_ret["return"]
        # print(json.dumps(d_bootstrapOccurrence, indent=4, sort_keys=True))

        d_ret = pipeline.innerLoop(
                    pipeline.bootstrap_thresholdDictionaryBuild,
                    threshold   = args.threshold,
                    log         = "Building threshold dictionaries...\n"
                )
        # print(json.dumps(pipeline._d_bootstrapThreshold, indent=4, sort_keys=True))

        d_ret = pipeline.innerLoop(
                    pipeline.bootstrap_filteredDictionaryBuild,
                    log         = "Building filtered dictionaries...\n"
                )
        # print(json.dumps(pipeline._d_bootstrapFiltered, indent=4, sort_keys=True))

        stage.exitCode(0)
        return True
    stage1.def_stage(f_stage1callback, obj=stage1, pipe=roigcibs)

    stage2 = stage.Stage(
                        name            = 'roigcibs-2-labelReader',
                        fatalConditions = True,
                        syslog          = True,
                        logTo           = '%s/roigcibs-2-labelReader.log'  % args.workingDir,
                        logTee          = True
                        )
    # stage4.def_preconditions(stage1.def_postconditions()[0], **stage1.def_postconditions()[1])
    stage2.def_preconditions(lambda **x: True)
    def f_stage2callback(**kwargs):
        for key, val in kwargs.iteritems():
            if key == 'obj':    stage       = val
            if key == 'pipe':   pipeline    = val

        os.chdir(pipeline._workingDir)

        d_ret = pipeline.innerLoop(
            pipeline.labelScript_process,
            log = "Writing and processing FreeSurfer tcl label script files...\n"
        )

        stage.exitCode(0)
        return True
    stage2.def_stage(f_stage2callback, obj=stage2, pipe=roigcibs)


    roigcibslog = roigcibs.log()
    roigcibslog('INIT: %s\n' % ' '.join(sys.argv))
    roigcibs.stage_add(stage0)
    roigcibs.stage_add(stage1)
    roigcibs.stage_add(stage2)
    # roigcibs.stage_add(stage3)
    # roigcibs.stage_add(stage4)
    roigcibs.initialize()

    roigcibs.run()
