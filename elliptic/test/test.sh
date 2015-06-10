#!/bin/sh

host=`hostname`

run_test() {
	for ((i = $1; i <= $2; i += $3))
	do
		echo `date` $i
		prog="elliptic $5"
		$prog -n $i >> $4
	done
}

#echo CG

#run_test 32 512 32 ${host}_p1_cg_s1.out "-c -p 1 -s 1"
#run_test 576 1024 64 ${host}_p1_cg_s1.out "-c -p 1 -s 1"

#run_test 32 512 32 ${host}_p1_cg_s4.out "-c -p 1 -s 4"
#run_test 576 1024 64 ${host}_p1_cg_s4.out "-c -p 1 -s 4"

#run_test 32 512 32 ${host}_p2_cg_s1.out "-c -p 2 -s 1"
#run_test 576 1024 64 ${host}_p2_cg_s1.out "-c -p 2 -s 1"

#run_test 32 512 32 ${host}_p2_cg_s2.out "-c -p 2 -s 2"
#run_test 576 1024 64 ${host}_p2_cg_s2.out "-c -p 2 -s 2"

#run_test 32 512 32 ${host}_p2_cg_s4.out "-c -p 2 -s 4"
#run_test 576 1024 64 ${host}_p2_cg_s4.out "-c -p 2 -s 4"

echo Jacobi
run_test 1024 1024 1 ${host}_p1_jac_s1.out "-j -p 1 -s 1"

echo GS
run_test 1024 1024 1 ${host}_p1_gs_s1.out "-g -p 1 -s 1"

echo GS
run_test 32 512 32 ${host}_p2_gs_s4.out "-g -p 2 -s 4"
run_test 1024 1024 1 ${host}_p1_gs_s2.out "-g -p 1 -s 2"
run_test 1024 1024 1 ${host}_p1_gs_s4.out "-g -p 1 -s 4"
run_test 1024 1024 1 ${host}_p2_gs_s1.out "-g -p 2 -s 1"
run_test 1024 1024 1 ${host}_p2_gs_s2.out "-g -p 2 -s 2"
run_test 1024 1024 1 ${host}_p2_gs_s4.out "-g -p 2 -s 4"

echo `date` Complete!
