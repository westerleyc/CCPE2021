# change number of threads
NUM_THREADS=8

export OMP_NUM_THREADS=$NUM_THREADS

set -e

g++ -std=c++11 -O3 -fopenmp sa_homog_VPR.cpp 

GRAPH=(
	'arf'
)

echo "benchmark, cost_mesh, cost_1-hop, time_mesh(s), time_1hop(s)" > results/results.csv 
for ((i = 0; i < ${#GRAPH[@]}; ++i)) do
    echo "list/"${GRAPH[i]}".in"
	./a.out "list/"${GRAPH[i]}".in" "names/"${GRAPH[i]}".in" ${GRAPH[i]} >> results/results.csv
done

rm a.out
