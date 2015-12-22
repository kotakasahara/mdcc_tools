#!/bin/bash
#$ -S /bin/bash
#$ -cwd
source /home/kotakasahara/.zshrc

/home/kotakasahara/local/kprml/kprml -mode test -feature 1 6 x -data-skip 1 -header-skip 1  -fn-data-table /home/kotakasahara/cal/transcript/ph08_ets1/analysis/analysis20131005/st04_3mfk/trj_cal08/nvt_fit_dt10_t.trr -format-data-table kktrajtrans
