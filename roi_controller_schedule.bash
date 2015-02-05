#!/bin/bash

# "include" the set of common script functions
source common.bash
declare -i Gi_verbose=1

# Column formatting (from common.bash)
G_LC=30
G_RC=50


G_SYNOPSIS="

  NAME

        roi_controller_schedule.bash

  SYNOPSIS

        roi_controller_schedule.bash                            \\
                                -r <subSampleList>              \\
                                -t <subSampleThreshold>         \\
                                -a <annotationList>             \\
                                -s <surfaceList>                \\
                                -h <hemiList>                   \\
                                -o <outputDirStem>              \\
                                -S <stages>                     \\
                                -g <groupList>                  \\
                                -c <curvFuncList>               \\
                                -C                              \\
                                -K                              \\
                                -R <computeNodeScratch>         \\
                                <SUBJ1> <SUBJ2> ... <SUBJn>

  DESC

        'roi_controller_schedule.bash' explicitly schedules each
        possible 'sub-run' of the roi analysis controller based on
        the pattern of input arguments.

  ARGS

    -r <subSampleList>
    A comma separated list of subSamples to create and analyze,
    eg. 1,2,3,4,5,6

    -t <subSampleThreshold>
    The threshold 0.0 ... 1.0 of resampling. 1.0 denotes a straight copy,
    0.5 denotes half of the original sample set.

    -a <annotationList>
    A comma separated list of annotations to analyze,
    eg. aparc.annot,aparc.a2009s.annot

    -s <surfaceList>
    A comma separeted list of surfaces to analyze,
    eg. smoothwm,pial

    -h <hemiList>
    A comma separeted list of hemispheres to analyze,
    eg. rh,lh

    -o <outputDirStem>
    The 'stem' name of the output directory, eg 4-dyslexia-run

    -S <stageList>
    A comma separated list of stages of the underlying roi_controller.py
    script to process. Note, this should start from stage 2.

    -g <groupList>
    The group of comparisons to process.

    -R <computeNodeScratch>
    The 'scratch' space directory in which this script will organize
    scheduler output (stdout and stderr) results from scheduled jobs.

    -C
    Clear existing ROItag directory. Since analysis run effectively
    'appends' to ROI hits, it is imperative to Clear existing directory
    trees from earlier runs, should they exist.

    -K
    Skip the intial resample step. This is necessary if an experiment is being
    run with an additional annotation file, for example.

    -N
    Skip remaining resample steps. This is necessary if part of the curvature
    analysis has already been completed and subsequent stages of roi_controller.py
    need to be run.

  HISTORY

  20-Aug-2014
  o Initial design and coding.

"

RESAMPLELIST="1,2,3"
RESAMPLETHRESHOLD="0.5"
ANNOTLIST="aparc.a2009s.annot,aparc.annot"
SURFACELIST="pial,smoothwm"
HEMILIST="lh,rh"
OUTSTEM="../4-exp-dyslexia-run"
STAGELIST="2,3,4"
GROUPLIST="13,26,24,45,46,56"
G_CNODESCRATCHPATH="~/scratch"
G_CURVFUNCLIST="H,K,K1,K2,C,BE,S,thickness"

let b_clearROI=0
let b_keepSamples=0
let b_skipRemainingRsample=0
G_LC=60
G_RC=30
ORIGARGS=$*

while getopts r:t:c:a:s:h:o:v:S:g:R:CKN option ; do
        case "$option"
        in
                r) RESAMPLELIST=$OPTARG         ;;
                t) RESAMPLETHRESHOLD=$OPTARG    ;;
                c) G_CURVFUNCLIST=$OPTARG       ;;
                a) ANNOTLIST=$OPTARG            ;;
                s) SURFACELIST=$OPTARG          ;;
                S) STAGELIST=$OPTARG            ;;
                h) HEMILIST=$OPTARG             ;;
                o) OUTSTEM=$OPTARG              ;;
                g) GROUPLIST=$OPTARG            ;;
                C) b_clearROI=1                 ;;
                K) b_keepSamples=1              ;;
                N) b_skipRemainingRsample=1     ;;
                R) G_CNODESCRATCHPATH=$OPTARG   ;;
                v) let Gi_verbose=$OPTARG       ;;
                \?) synopsis_show
                    exit 0;;
        esac
done

verbosity_check
topDir=$(pwd)


# Remove commas from STAGELIST
STAGELIST=$(echo $STAGELIST | sed 's/\,//g')

shift $(($OPTIND - 1))
SUBJLIST=$*
b_SUBJLIST=$(echo $SUBJLIST | wc -w)

stage_stamp "${G_SELF} $ORIGARGS" "${G_SELF}.log" "INIT"

if (( b_clearROI )) ; then
    dirHere=$(pwd)
    for RESAMPLE in $(echo $RESAMPLELIST | tr ',' ' '); do
        for ANNOT in $(echo $ANNOTLIST | tr ',' ' '); do
            if [[ -d  ${OUTSTEM}${RESAMPLE}/groupCurvAnalysis/$ANNOT ]] ; then
                lprint "Deleting ROItag-$G_CURVFUNCLIST tree in sample $RESAMPLE, $ANNOT..."
                cd ${OUTSTEM}${RESAMPLE}/groupCurvAnalysis/$ANNOT
                echo "rm -fr ROItag-$G_CURVFUNCLIST" | sh
                rprint "[ ok ]"
                cd $dirHere
            fi
        done
    done
fi

# First do the resampling/setup... there's a clash when run in parallel
# so need to run *once* first per RESAMPLE with stage 0, and then all
# the annotations for stages 1 and 2.
if (( !b_keepSamples )) ; then
    for RESAMPLE in $(echo $RESAMPLELIST | tr ',' ' '); do
        CMD="/neuro/users/rudolphpienaar/src/devel/roi_tag/rsample.py\
            -o ${OUTSTEM}${RESAMPLE} -r $RESAMPLETHRESHOLD -s 0 -a $ANNOT $SUBJLIST
            "
        echo "$CMD"
        eval "$CMD"
    done
fi

if (( ! b_skipRemainingRsample )) ; then
    for RESAMPLE in $(echo $RESAMPLELIST | tr ',' ' '); do
        for ANNOT in $(echo $ANNOTLIST | tr ',' ' '); do
            CMD="/neuro/users/rudolphpienaar/src/devel/roi_tag/rsample.py\
                -o ${OUTSTEM}${RESAMPLE} -r $RESAMPLETHRESHOLD -s 12 -a $ANNOT $SUBJLIST
                "
            echo "$CMD"
            eval "$CMD"
        done
    done
fi

# Now the remaining curvature analysis stages can be each be run in
# parallel, skipping stages '01' of roi_controll.py since the resampling
# has already been done.
#
# Each stage spawns a set of child processes, and the crun object will
# block until all children have finished before looping to the next
# stage.
#
# NB! When running through a scheduler, pay particular attention to the
# output and working directory concepts!
THISDIR=$(pwd)

CMD=$(echo "                        \
/neuro/users/rudolphpienaar/src/devel/roi_tag/roi_controller.py                       \
    -w $THISDIR                     \
    -r $RESAMPLELIST                \
    -t $RESAMPLETHRESHOLD           \
    -v 1                            \
    -a $ANNOTLIST                   \
    -o $OUTSTEM                     \
    -s $STAGELIST                   \
    --curvFunc $G_CURVFUNCLIST      \
    --pval le1                      \
    --surface $SURFACELIST          \
    --hemi $HEMILIST                \
    --statFunc ptile-raw            \
    --group $GROUPLIST              \
    --schedulerStdOutDir $G_CNODESCRATCHPATH/$$-roi_controller.jobout    \
    --schedulerStdErrDir $G_CNODESCRATCHPATH/$$-roi_controller.joberr    \
    $SUBJLIST " | tr '\n' ' ' | sed 's/ \+/ /g')
echo "$CMD"
eval "$CMD"

shut_down 0
