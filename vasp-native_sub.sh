#!/bin/bash

EMAIL=AUTO

# You shouldn't need to change anything else. Run this script with -h for usage.
#-------------------------------------------------------------------------------

if [ "$EMAIL" = AUTO ]; then
  EMAIL="$(whoami)@vt.edu"
fi

QUEUE=normal_q
LIST=newriver
NODES=1
CPUSPER=24
WALLTIME="48:00:00"
SECS=600

while getopts  "c:ad:e:hj:mn:pq:sw:" flag
do
  case $flag in
    h) echo "`basename $0` Version $VERSION for PBS (torque) Copyright 2005 - Peter Larkowski"
       echo
       echo "Usage: `basename $0` [-h] [-c CPUSPER] [-j JOBNAME ] [-m [-d SEC]]"
       echo "                     [-n NUM_NODES ] [-s] [-w WALLTIME ] JOB_FILE(S)"
       echo "   -h - Print this message"
       echo "   -c - CPUSPER is the number of cpus per node the job will use (default 2)"
       echo "   -j - Set Job name to JOBNAME (default is JOB_FILE)"
       echo "   -m - Turn on memory usage logging."
       echo "   -d - Take a memory reading every SEC seconds.  Default is every 10 min"
       echo "        (600 seconds).  Ignored if -m is not specified"
       echo "   -n - NUM_NODES is the number of nodes to run on (default is 1)"
       echo "   -p - Turns on the preempt option"
       echo "   -q - QUEUE is the queue to submit to (default is normal_q)"
       echo "   -l - LIST is the grouplist to submit to (default is blueridge)"
       echo "   -s - Submit job after generating the script"
       echo "   -w - WALLTIME is the walltime (default 100:00:00)" 
       echo "   -a - increase memory"
       exit 0
     ;;
    j) JOBNAME=$OPTARG
     ;;
    d) if [ "$(echo $OPTARG | grep -v [^0-9])" = "" ]; then
         echo "`basename $0`: -d requires an integer argument"
         exit 1
       else
         SECS=$OPTARG
       fi
     ;;
    n) if [ "$(echo $OPTARG | grep -v [^0-9])" = "" ]; then
         echo "`basename $0`: -n requires an integer argument"
         exit 1
       else
         NODES=$OPTARG
       fi
     ;;
    c) if [ "$(echo $OPTARG | grep -v [^0-9])" = "" ]; then
         echo "`basename $0`: -c requires an integer argument"
         exit 1
       else
         CPUSPER=$OPTARG
       fi
     ;;
    m) MEMORY=yes
     ;;
    p) PAY=yes
     ;;
    s) SUBMIT=yes
     ;;
    w) if [ "$(echo $OPTARG | grep -v [^0-9:])" = "" ]; then
         echo "`basename $0`: -w requires an argument that looks like: XXX:XX:XX"
         exit 1
       else
         WALLTIME=$OPTARG
       fi
     ;;
    e) EMAIL=$OPTARG
     ;;
    ?) echo "`basename $0`: Unknown option $flag try -h"
       exit 1
     ;;
  esac
done
shift $(($OPTIND - 1))

PWD=`pwd`

if [ "$*" = "" ]; then
  echo "You must give me at least 1 job file to work with..."
  echo "`basename $0` -h for more info"
  exit 1
fi

if [ "$MEMORY" = yes ]; then
    MEMORY=":highmem"
else
    MEMORY=""
fi
if [ "$PAY" = yes ]; then
    QUEUE="hxin_lab"
    LIST="hxin_lab"
fi

for i in $*; do
  BASE="$(basename $i .drun)"

  if [ "$JOBNAME" = "" ]; then
    JOBNAME=$BASE
  fi
  echo $JOBNAME "-" $i
 
  cat << EOF > $BASE.qsub
#!/bin/bash
#PBS -S /bin/bash
#PBS -N $JOBNAME
#PBS -l nodes=$NODES:ppn=$CPUSPER$MEMORY,walltime=$WALLTIME
#PBS -q $QUEUE
#PBS -M $EMAIL
#PBS -W group_list=$LIST
#PBS -m ae
#PBS -o std.o
#PBS -e std.e
#PBS -V
#PBS -A DFT_XinLab

cd $PWD

export MACHINE_COUNT=`echo $NODES*$CPUSPER | bc`

EOF

  if [ "$DEBUGM" = yes ]; then
    cat << EOF >> $BASE.qsub
JOBID="\$(echo \$PBS_JOBID | cut -f1-2 -d.)"

report-mem-usage() {
  while true; do
    date >> memusage-$JOBNAME.log
    for j in \`cat \$PBS_NODEFILE | uniq\`; do
      echo \$j >> memusage-$JOBNAME.log
      ssh \$j free >> memusage-$JOBNAME.log
    done
    sleep $SECS
    if [ "\$(qstat | grep \$JOBID)" = "" ]; then
      exit 0
    fi
  done
}

report-mem-usage &

EOF
  fi

  cat << EOF >> $BASE.qsub
mpirun -np \$MACHINE_COUNT -machinefile \$PBS_NODEFILE /bin/bash -c "ulimit -s unlimited; vasp_std"

EOF

chmod 755 $BASE.qsub

if [ "$SUBMIT" = yes ]; then
  qsub $BASE.qsub
else
  echo "Now submit the job with \"qsub $BASE.qsub\""
fi

done 
