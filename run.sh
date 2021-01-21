# change number of threads
NUM_THREADS=8

export OMP_NUM_THREADS=$NUM_THREADS

set -e

g++ -std=c++11 -O3 -fopenmp sa_homog_VPR.cpp 

GRAPH=(
	'mac'
	'simple'
	'horner_bs'
	'mults1'
	'arf' 
	'conv3' 
	'motion_vec'
	'fir2'
	'fir1'
	'fdback_pts'
	'k4n4op'
	'h2v2_smo'
	'cosine1' 
	'ewf' 
	'Cplx8' 
	'Fir16'
	'cosine2' 
	'FilterRGB'
	'collapse_pyr'
	'interpolate'
	'w_bmp_head'
	'matmul'
	'invert_matrix'
)

echo "benchmark, nodes, edges, gridSize, cost_1-hop, time_1hop(s), swapCount" > results/results.csv 
for ((i = 0; i < ${#GRAPH[@]}; ++i)) do
    echo "list/"${GRAPH[i]}".in"
	./a.out "list/"${GRAPH[i]}".in" "names/"${GRAPH[i]}".in" ${GRAPH[i]} >> results/results.csv
done

rm a.out
