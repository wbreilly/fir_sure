#!/bin/bash -l 
#
#SBATCH --job-name=fir_10_30_18
#SBATCH --output=/home/wbreilly/sms_scan_crick/cluster_fir_data_10_27_18/logs/fir_10_30_18.%j.%N.out
#SBATCH --error=/home/wbreilly/sms_scan_crick/cluster_fir_data_10_27_18/logs/fir_10_30_18.%j.%N.err
#SBATCH --nodes=1
#SBATCH -c 34
#SBATCH -p bigmemm
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=3000 # 3gb % spm max mem set at 2.1gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wbreilly@ucdavis.edu

date 
hostname
module load matlab/R2016a # load matlab module 

matlab -nodisplay < parbatch_fir.m # run matlab script 
