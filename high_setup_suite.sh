#!/bin/bash


RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color





##################################################
  #   Run from monc main directory   # 
##################################################
dbase=$(pwd)
cd $dbase
echo ""
echo -e "   Setting up suite for ${GREEN}${dbase}${NC}..."
echo ""
##################################################




##################################################
   #   Options
##################################################

# Essential:
project='n02-REVCON'

# Requires: 
mcf_template=config_template_weak.mcf
submit_template=submonc_template.sb

stag=cdum       # Run type prefix label
# Default: 
ncols=32                # number of model columns per MONC processor, in x and in y
a2cores_per_node=128    # number of processors per node on ARCHER2


# Processor (resource size)  to set up (2**19 = 524288 cores, 4096/5276 ARCHER2 nodes)
# range of 2**exponent (2**7=128 --> one ARCHER2 node)
min_pe=0
max_pe=14

# Set up number of columns per MONC (problem size), in x and in y
# range of 2**exponent (2**5=32 --> default choice?)
#   - each step is 4x horizontal domain on a MONC
#   - NOTE: MONC will fail immediately on a fully utilised node for the 400m RCE case
#           with mpio=15 and ncols > ~100.
#           It will run without any IO sampling for ncols ~= 90.
min_ne=4
max_ne=4


# number of hours to run monc per submission cycle
internal=1

# number of timsteps to complete
#  - A value of -1 will set the model to run to its termination time 
set_nn=10000

# set up number of columns per MONC
for nnc in $(seq $min_ne $max_ne) ; do

  ncols=$((2**${nnc}))

  ntag=$(printf "%07d" ${ncols})

# Loop over core count exponent
for pnc in $(seq $min_pe $max_pe) ; do

  ncores=$((2**${pnc}))

  # Determine moncs_per_io_server
  # For small jobs
  if [[ ${ncores} -lt 16 ]] ; then
    mpio=$ncores
    ncores=$((${ncores} + 1))
    nmoncs=${mpio}
  else
    mpio=15
    nmoncs=$((${ncores} - ${ncores} / 16))
  fi

  # Configure number of cores label
  ctag=$(printf "%07d" ${ncores})


  # Determine node count
  if [[ ${ncores} -lt ${a2cores_per_node} ]] ; then
    nselect=1
    ntasks=${ncores}
  else
    nselect=$((${ncores} / ${a2cores_per_node}))
    ntasks=${a2cores_per_node}
  fi



  # get largest factor pair (nearly square domain)
  bbb=($(python factors.py ${nmoncs}))
  xmoncs=${bbb[0]}
  ymoncs=${bbb[1]}
  xsize=$(( ${xmoncs} * ${ncols} ))
  ysize=$(( ${ymoncs} * ${ncols} ))

#  case_id=${stag}_${ctag}
  case_id=${stag}_c${ntag}_p${ctag}




  # Process more jobs at once under taskfarm
  if [[ ${nselect} -le 16 ]] ; then
    queue='taskfarm'
  else
    queue='standard'
  fi




  echo $ncores $nselect $case_id $mpio $nmoncs $xmoncs $ymoncs $xsize $ysize

  # Save some parameters
  qname=${case_id}
  cdir=$dbase/suite/$case_id

  # Set up directories and support files
  if [ ! -d $cdir ]; then
    mkdir $cdir
  else
    echo ""
    echo -e "[${RED}WARN${NC}]---------------------------------------------------------------"
    echo -e "[${YELLOW}WARN${NC}] The directory for ${GREEN}${case_id}${NC} already exists."
    echo -e "[${RED}WARN${NC}] Do you wish to overwrite, deleting the entire case"
    echo -e "[${YELLOW}WARN${NC}] and any data within that directory?"
    echo -e "[${RED}WARN${NC}]"
    echo -e "[${YELLOW}WARN${NC}]  - ${BLUE}Enter the number corresponding to your response.${NC} - "
    echo -e "[${RED}WARN${NC}]---------------------------------------------------------------"
    select yn in "Yes" "No"; do
      case $yn in
        Yes ) echo 'You selected Yes.  Removing...'; rm -rf $cdir/*; break;;
        No ) echo 'You selected No.  Exiting...'; exit 1 ;;
      esac
    done
  fi

  echo ""
  echo Creating $cdir
#  echo   Resolution: $res
  echo   NGPTS:      $xsize, $ysize
  echo   Nodes:      $nselect
  echo   Walltime:   $internal
  echo   Job Name:   $qname
  echo ""

  mkdir $cdir/diagnostic_files
  diagdir=suite/$case_id/diagnostic_files
  mkdir $cdir/checkpoint_files
  ckptdir=suite/$case_id/checkpoint_files
  mkdir $cdir/monc_stdout
  stddir=suite/$case_id/monc_stdout
  vncp=suite/$case_id/version
  cp $cdir/../${submit_template} $cdir/csubmonc.sb
  pbs=suite/$case_id/csubmonc.sb
  cp $cdir/../${mcf_template} $cdir/config.mcf
  config=suite/$case_id/config.mcf
  cp $dbase/misc/continuation.sh $cdir/continuation.sh
  conscript=suite/$case_id 

  cp $dbase/version $cdir/version

  # Link in the start dump file
#  startdump=`/bin/ls -rt1 $dnc/checkpoint_files | tail -1`
#  if [ -z "$startdump" ]; then
#    echo ""
#    echo -e "${RED} --WARNING--${NC}"
#    echo "  There were no checkpoint files to draw from.  Go find one."
#    echo -e "${RED} --WARNING--${NC}"
#  else
#    echo ""
#    echo "Linking $dnc/checkpoint_files/$startdump to $dbase/$ckptdir/${case_id}_dump_startdump.nc"
#
#    ln -s $dnc/checkpoint_files/$startdump $dbase/$ckptdir/${case_id}_dump_startdump.nc
#  fi

  # Swap items in files
  # Submission script pbs:
  # #PBS -N add_qname
    sed -i -e 's/add_qname/'"$qname"'/g' $pbs

  # #PBS -l select=add_select
    sed -i -e 's/add_select/'"$nselect"'/g' $pbs
    sed -i -e 's/add_np/'"$ncores"'/g' $pbs

  # #PBS -l walltime=add_walltime
    sed -i -e 's/add_walltime/'"$internal:40:00"'/g' $pbs

  # #PBS -o add_odir
    sed -i -e 's|add_odir|suite/'"$case_id"'/|g' $pbs

  # #PBS -P add_project
    sed -i -e 's/add_project/'"$project"'/g' $pbs

  # #PBS -q add_queue
    sed -i -e 's/add_queue/'"$queue"'/g' $pbs

  # export SUBMISSION_SCRIPT_NAME=add_scriptname
    sed -i -e 's|add_scriptname|'"$pbs"'|g' $pbs

  # export TESTCASE=add_mcf
    sed -i -e 's|add_mcf|'"$config"'|g' $pbs

  # export STDOUT_DIR=add_stdout_dirname
    sed -i -e 's|add_stdout_dirname|'"$stddir"'|g' $pbs

  # cp ./version add_vncp
    sed -i -e 's|add_vncp|'"$vncp"'|g' $pbs

  # export CP_DIR=add_checkpoint_dirname
    sed -i -e 's|add_checkpoint_dirname|'"$ckptdir"'|g' $pbs

  # export RUN_NAME=add_runname
    sed -i -e 's|add_runname|'"$case_id"'_dump_|g' $pbs

  # . add_condir/continuation.sh
    sed -i -e 's|add_condir|'"$conscript"'|g' $pbs


  # Configuration script items:
  # checkpoint_file="add_checkpoint_dir/add_run_name.nc"     # Checkpoint file location and prefix
    sed -i -e 's|add_checkpoint_dir|'"$ckptdir"'|g' $config
    sed -i -e 's|add_run_name|'"$case_id"'_dump|g' $config

  # walltime_limit=add_internal_walltime          # Internal wall clock time limit on simulation [hh:mm:ss]
    sed -i -e 's|add_internal_walltime|'"$internal:00:00"'|g' $config

  # x_size=add_xsize
  # y_size=add_ysize
    sed -i -e 's|add_xsize|'"$xsize"'|g' $config
    sed -i -e 's|add_ysize|'"$ysize"'|g' $config

  # moncs_per_io_server=add_mpio
    sed -i -e 's|add_mpio|'"$mpio"'|g' $config

  # nn_timesteps=add_nn
    sed -i -e 's|add_nn|'"$set_nn"'|g' $config

  # z_size=add_z_size
#    sed -i -e 's|add_nz1|'"${nz1[cnc]}"'|g' $config
#    sed -i -e 's|add_nz2|'"${nz2[cnc]}"'|g' $config
#    sed -i -e 's|add_z_size|'"${nz3[cnc]}"'|g' $config

  # dxx=add_res
  # dyy=add_res
#    sed -i -e 's|add_res|'"$res"'|g' $config

  # f_force_pl_theta=add_fcg, add_fcg, 0.0, 0.0
#    sed -i -e 's|add_fcg|-'"$fcg_id"'|g' $config

  # diagnostic_file="add_diagnostic_file"
  #     diagnostic_files/case_casim_dg.nc
    sed -i -e 's|add_diagnostic_file_1file|'"$diagdir/${case_id}_1file.nc"'|g' $config


  echo "Done."
  echo ""

done # pnc loop over core count exponent

done # cnc loop over compiler tags

exit 0


