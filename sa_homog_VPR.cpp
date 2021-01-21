#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <ctime>
#include <chrono>
#include <algorithm> 
#include <fstream>
#include <omp.h>

#define NGRIDS 10
#define EMPTY 255

using namespace std;
using namespace std::chrono;

long long int results[NGRIDS][4];
double randomvec[1000000];

//Tables with costs for mesh and 1hop
vector<vector<int>> tablemesh;
vector<vector<int>> table1hop;

/* GLOBAL VARIABLES INIT */
    int bestCost;
    int idxMinCost;

    //graph information
    int nodes, edges, nIO;
    int *h_edgeA, *h_edgeB;
    vector<int> A;
    int *v, *v_i;
    bool *io;
    map<int,string> names;

    //placement information
    int dim, gridSize;
    int **grid;             
    int cost = 100000;
    int **positions;
/* GLOBAL VARIABLES END */

void readInput(char* arqname){
    FILE *fptr;
	fptr = fopen(arqname,"r");
    if (fptr == NULL) {
        printf("Error! opening file\n");
        exit(1);
    }

    //INPUT READ
    int c = 0, n1, n2;
    //NUMERO DE GRIDS
    int nGrids = NGRIDS; 
    if(fscanf(fptr, "%d %d", &n1, &n2)){
        nodes = n1;
        edges = n2;
    }
    h_edgeA = new int[edges];
    h_edgeB = new int[edges];
    v = new int[nodes];
    v_i = new int[nodes];
    io = new bool[nodes];
    for(int i=0; i<nodes; i++){
        v[i]=0; v_i[i]=0;
        io[i] = false;
    }
    c=edges;
    int counter=0;  
    while(c>0) {
        if(fscanf(fptr, "%d %d", &n1, &n2)){
            h_edgeA[c-1] = n1;
            h_edgeB[c-1] = n2;
        }
        v[n1]++;
        if(n1!=n2) v[n2]++;
        c--;
    }
    while(fscanf(fptr, "%d", &n1) != EOF){
        io[n1] = true;
        counter++;
    }
    nIO = counter;

    //Complete graph information
    for(int i=1; i<nodes; i++){
        v_i[i] = v_i[i-1] + v[i-1];
    }
    for(int i=0; i<nodes; i++){
        for(int j=0; j<edges; j++){
            if (h_edgeA[j] != h_edgeB[j]) {
                if(h_edgeA[j]==i) A.push_back(h_edgeB[j]);
                if(h_edgeB[j]==i) A.push_back(h_edgeA[j]);
            } else {
                if(h_edgeA[j]==i) A.push_back(h_edgeB[j]);
            }
        }
    } 
}

void readNames(string arqname){
    ifstream fptr;
	fptr.open(arqname);
    string n1, n2;
    while(!fptr.eof()){
        fptr >> n1;
        fptr >> n2;
        int n3 = stoi(n1);
        names[n3] = n2;
    }
}

void buildMatrices(int dim){
    table1hop.resize(dim);
    tablemesh.resize(dim);
    for(int i=0; i<dim; i++){
        table1hop[i].resize(dim);
        tablemesh[i].resize(dim);
    }

    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            int jfrom = 0;
            int ifrom = 0;
            int distX = abs(jfrom-j);
            int distY = abs(ifrom-i); 
            tablemesh[i][j] = max(1,distY+distX);
            table1hop[i][j] = max(1,distX/2+distX%2+distY/2+distY%2);
        }
    }
}

void printGraphInfo(){
    printf("---------------GRAPH INFORMATION---------------\n");
    printf("NODES: %d EDGES: %d\n", nodes, edges);
    for (int i = 0; i < edges; ++i) printf("%d -> %d\n", h_edgeA[i], h_edgeB[i]);    
    printf("v: ");
    for (int i = 0; i < nodes; ++i) printf("%d ", v[i]);
    printf("\n");
    printf("v_i: ");
    for (int i = 0; i < nodes; ++i) printf("%d ", v_i[i]);
    printf("\n");
    printf("A: ");
    for (int i = 0; i < A.size(); ++i) printf("%d ", A[i]);
    printf("\n");
    printf("ios: ");
    for (int i = 0; i < nodes; ++i) if(io[i]) printf("%d ", i);
    printf("\n");
    printf("-----------------------------------------------\n\n");
}

void printGrid(int idx){
    printf("grid:\n");
    for (int i = 0; i < gridSize; ++i){
        if(i%dim==dim-1) printf("%d\n", grid[idx][i]);
        else printf("%d\t", grid[idx][i]);
    }
    printf("\n");
}

int gridCost(int iresult, int arch){
    int cost_ = 0, increment=0, distManhattanI, distManhattanJ;
    for(int k=0; k<edges; k++){
        int ifrom = positions[iresult][h_edgeA[k]]/dim; 
        int jfrom = positions[iresult][h_edgeA[k]]%dim; 
        int ito = positions[iresult][h_edgeB[k]]/dim; 
        int jto = positions[iresult][h_edgeB[k]]%dim;
        
        distManhattanJ = abs(jto - jfrom);
        distManhattanI = abs(ito - ifrom);

        if (arch == 0)
            increment = tablemesh[distManhattanI][distManhattanJ];
        else if (arch == 1) {
            increment = table1hop[distManhattanI][distManhattanJ];                                 
        } 
        cost_ += increment;
    }
    return cost_;
}

void annealing2(int iresult, int arch, bool* borders, bool* corners){
    int *localGrid, *localPositions;
    int currentCost = cost, nextCost, increment=0, distManhattanI, distManhattanJ;
    long long int swapCount = 0;
    localGrid = new int[gridSize];
    for(int i=0; i<gridSize; i++) localGrid[i] = grid[iresult][i];
    localPositions = new int[nodes];
    for(int i=0; i<nodes; i++) localPositions[i] = positions[iresult][i];    
    //random vector index
    double random, valor;
    int randomctrl = 0;
    double T=100;
    auto start = high_resolution_clock::now();
    while(T>=0.00001){
        for(int i=0; i<gridSize; i++){
            for(int j=i+1; j<gridSize; j++){                
                //if we're looking at 2 empty spaces, skip                   
                if(localGrid[i]==EMPTY && localGrid[j]==EMPTY)
                    continue;
                if(corners[i]==true || corners[j]==true) 
                    continue;
                //swap only if both nodes are on the border or both nodes are internal
                if(borders[i]==true && borders[j]==false)
                    continue;
                if(borders[i]==false && borders[j]==true)
                    continue;

                int node1 = localGrid[i], node2 = localGrid[j];
                nextCost = currentCost;
                int old1, old2;
                bool test = false;

                if(true){
                    //remove cost from object edges                
                    if(node1!=EMPTY){
                        for(int i=0; i<v[node1]; i++){
                            //nextCost -= calcIncrement(node1);
                            int ifrom = localPositions[node1]/dim; 
                            int jfrom = localPositions[node1]%dim; 
                            int ito = localPositions[A[v_i[node1]+i]]/dim; 
                            int jto = localPositions[A[v_i[node1]+i]]%dim;
                            distManhattanJ = abs(jto - jfrom);
                            distManhattanI = abs(ito - ifrom);
                            if (arch == 0)
                                increment = tablemesh[distManhattanI][distManhattanJ];
                            else if (arch == 1) {
                                increment = table1hop[distManhattanI][distManhattanJ];     
                            } 
                            nextCost -= increment; 
                        }
                    }
                    if(node2!=EMPTY){
                        for(int i=0; i<v[node2]; i++){
                            int ifrom = localPositions[node2]/dim; 
                            int jfrom = localPositions[node2]%dim; 
                            int ito = localPositions[A[v_i[node2]+i]]/dim; 
                            int jto = localPositions[A[v_i[node2]+i]]%dim; 
                            distManhattanJ = abs(jto - jfrom);
                            distManhattanI = abs(ito - ifrom);

                            if (arch == 0)
                                increment = tablemesh[distManhattanI][distManhattanJ];
                            else if (arch == 1) {
                                increment = table1hop[distManhattanI][distManhattanJ];     
                            }
                            nextCost -= increment; 
                        }
                    }
                    //swap positions
                    old1 = i;
                    old2 = j;
                    if(node1!=EMPTY) localPositions[node1] = old2;
                    if(node2!=EMPTY) localPositions[node2] = old1;
                    localGrid[j] = node1;
                    localGrid[i] = node2;
                    //recalculate cost
                    if(node1!=EMPTY){
                        for(int i=0; i<v[node1]; i++){
                            int ifrom = localPositions[node1]/dim; 
                            int jfrom = localPositions[node1]%dim; 
                            int ito = localPositions[A[v_i[node1]+i]]/dim; 
                            int jto = localPositions[A[v_i[node1]+i]]%dim; 
                            distManhattanJ = abs(jto - jfrom);
                            distManhattanI = abs(ito - ifrom);

                            if (arch == 0)
                                increment = tablemesh[distManhattanI][distManhattanJ];
                            else if (arch == 1) {
                                increment = table1hop[distManhattanI][distManhattanJ];     
                            }
                            nextCost += increment;                         
                        }
                    }
                    if(node2!=EMPTY){
                        for(int i=0; i<v[node2]; i++){
                            int ifrom = localPositions[node2]/dim; 
                            int jfrom = localPositions[node2]%dim; 
                            int ito = localPositions[A[v_i[node2]+i]]/dim; 
                            int jto = localPositions[A[v_i[node2]+i]]%dim; 
                            distManhattanJ = abs(jto - jfrom);
                            distManhattanI = abs(ito - ifrom);

                            if (arch == 0)
                                increment = tablemesh[distManhattanI][distManhattanJ];
                            else if (arch == 1) {
                                increment = table1hop[distManhattanI][distManhattanJ];     
                            }
                            nextCost += increment;   
                        }
                    }
                }
                //parameter for annealing probability
                valor = exp(-1*(nextCost - currentCost)/T);
                //random number between 0 and 1
                random = randomvec[randomctrl];
                randomctrl++;
                if(randomctrl==1000000) randomctrl=0;

                //if cost after changes is less than before or if cost is higher but we're in the annealing probability range, return
                if(nextCost <= currentCost || random <= valor){
                    currentCost = nextCost;
                    swapCount++;
                }
                //else, undo changes and stay with previous cost
                else{
                    if(node1!=EMPTY) localPositions[node1] = old1;
                    if(node2!=EMPTY) localPositions[node2] = old2;
                    localGrid[j] = node2;
                    localGrid[i] = node1;
                }
            }
            if(cost==edges) break;             
            T*=0.999;
        }   
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    for(int i=0; i<gridSize; i++) grid[iresult][i] = localGrid[i];
    for(int i=0; i<nodes; i++) positions[iresult][i] = localPositions[i]; 
    
    results[iresult][1] = currentCost;
    results[iresult][2] = duration.count()/1000;
    results[iresult][3] = swapCount;

    if(currentCost<=bestCost){
        bestCost = currentCost;
        idxMinCost = iresult;
    }
    delete localPositions;
    delete localGrid;
}

int main(int argc, char** argv) {
    srand (time(NULL));
    string arqName(argv[1]);
    string arqName2(argv[2]);
    for(int i=0;i<1000000;i++){
        randomvec[i] = (double)rand() / (double)(RAND_MAX);
    }

    //Read info from files
    readInput(argv[1]);
    //readNames(arqName2);

    //Precalculations
    dim = ceil(sqrt(nodes-nIO))+2;
    int internals = nodes-nIO;
    int freeBorders = 4*dim - 8;
    int freeIns = (dim-2)*(dim-2);
    while(nIO > freeBorders || internals > freeIns){
        dim++;
        freeBorders = 4*dim - 8;
        freeIns = (dim-2)*(dim-2);
    } 

    buildMatrices(dim);
    gridSize = dim*dim;

    //Results setup for mesh and hop architectures
    const int num_arch = 2;
    double time_total[num_arch] = {0.0, 0.0};
    int cost_min[num_arch] = {-1,-1};
    long int numberSwaps[num_arch] = {0,0};

    //Memory allocation for the grids
    grid = new int*[NGRIDS];
    for(int i=0; i<NGRIDS; i++) grid[i] = new int[gridSize];
    positions = new int*[NGRIDS];
    for(int i=0; i<NGRIDS; i++) positions[i] = new int[nodes];

    //Pattern for io placement on borders
    bool *borders, *corners;   
    borders = new bool[gridSize];
    corners = new bool[gridSize];
    for(int i=0; i<gridSize; i++){
        borders[i] = false;
        corners[i] = false;
        if(i/dim==0) borders[i] = true;
        else if(i/dim==dim-1) borders[i] = true;
        if(i%dim==0) borders[i] = true;
        else if(i%dim==dim-1) borders[i] = true;

        if(i==0 || i==dim-1 || i==gridSize-dim || i==gridSize-1) corners[i]=true;
    }

    //Output results
    ofstream ofile, ofile1;
    char buffer[strlen(argv[1])+20] = "\0";
    char beg[15] = "results/mesh/"; 
    char out[7] = ".out";
    strcat(buffer,beg);
    strcat(buffer,argv[3]);
    strcat(buffer,out);        

    char buffer1[strlen(argv[1])+20] = "\0";
    char beg1[15] = "results/1hop/"; 
    strcat(buffer1,beg1);
    strcat(buffer1,argv[3]);
    strcat(buffer1,out);
    ofile1.open(buffer1);

    char buffer2[strlen(argv[1])+20] = "\0";
    char beg2[15] = "data/net/dac/";
    char out2[7] = ".net"; 
    strcat(buffer2,beg2);
    strcat(buffer2,argv[3]);
    strcat(buffer2,out2);
    string arqnameX(buffer2);

    //k stands for type of architecture (mesh or 1hop)
    for (int k = 1; k < num_arch; ++k) {
        //Fill the grids
        for(int i=0; i<NGRIDS; i++){
            for(int j=0; j<gridSize; j++){
                grid[i][j] = EMPTY;
            }
        }
        for(int i=0; i<NGRIDS; i++){
            bool *setNodes;
            setNodes = new bool[nodes];
            for(int j=0; j<nodes; j++){
                setNodes[j] = false;
            }
            //Place IO nodes first
            for(int j=0; j<nodes; j++){
                if(io[j]==false) continue;
                for(int t=0; t<gridSize; t++){
                    if(borders[t]==true && setNodes[j]==false && grid[i][t]==EMPTY && corners[t] == false){
                        grid[i][t]=j;
                        setNodes[j] = true;
                        break;
                    }
                }
            }
            //Place other nodes            
            for(int j=0; j<nodes; j++){
                if(io[j]==true || setNodes[j]==true) continue;
                for(int t=0; t<gridSize; t++){
                    if(grid[i][t]==EMPTY && borders[t]==false){
                        grid[i][t]=j;
                        break;
                    }
                }
            }
            //Shuffle
            for(int j=0; j<gridSize; j++){
                for(int t=0; t<gridSize; t++){
                    double aleat = ((rand()%1000)/1000.0);
                    if(aleat > 0.2 && (borders[j]==true && borders[t]==true && corners[t] == false && corners[j] == false)){
                        int temp = grid[i][j];
                        grid[i][j] = grid[i][t];
                        grid[i][t] = temp;
                    }
                    if(aleat > 0.2 && (borders[j]==false && borders[t]==false)){
                        int temp = grid[i][j];
                        grid[i][j] = grid[i][t];
                        grid[i][t] = temp;
                    }  
                }
            }
        }

        //Fill positions vector
        for(int t=0; t<NGRIDS; t++){
            for(int i=0; i<gridSize; i++){
                for(int j=0; j<nodes; j++){
                    if(j==grid[t][i]) {
                        positions[t][j]=i;
                    }
                }
            }
        }

        //Annealing Step
        auto start = high_resolution_clock::now();
        #pragma omp parallel for
        for (int i = 0; i < NGRIDS; ++i) {
            cost = gridCost(i, k);
            results[i][0] = cost;
            if(cost == 0){
                results[i][1] = cost;
                results[i][2] = 0;
                results[i][3] = 0; 
                continue;
            }
            annealing2(i, k, borders, corners);
            if (results[i][3] > numberSwaps[k]) numberSwaps[k] = results[i][3];
        }
        auto stop = high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = (stop-start);
        time_total[k] = duration.count()/1000;

        vector<int> costs;    
        int totalTime=0;
        int cost_=0;
        for(int i=0; i<NGRIDS; i++){
            cost_ = gridCost(i,k);
            costs.push_back(cost_);
            // if(results[i][1]>=0) costs.push_back(results[i][1]);
            totalTime += results[i][2];
        }
        cost_min[k] = *min_element(costs.begin(), costs.end());

        //printing results
        for(int i=0; i<NGRIDS; i++){
            //ofile1 << gridCost(i,k) << endl;
            //ofile1 << results[i][1] << endl;
            for(int j=0; j<gridSize; j++){
                if(j == gridSize-1) ofile1 << grid[i][j] << endl;
                else ofile1 << grid[i][j] << " ";
            }
        }

        //DEBUG (GRIDS)
        // for(int i=0; i<NGRIDS; i++){
        //     ofile1 << "gridCost: "  << gridCost(i,k) << "results: " << results[i][1] << endl;
        //     for(int j=0; j<gridSize; j++){
        //         // if(j == gridSize-1) ofile1 << grid[i][j] << endl;
        //         // else ofile1 << grid[i][j] << " ";
        //         if(j%dim == dim-1) ofile1 << grid[i][j] << endl;
        //         else ofile1 << grid[i][j] << "\t";
        //     }
        //     ofile1 << endl;
        // }     


        // if(k==0){
        //     // printGrid(0);
        //     ofile.open(buffer);            
        //     ofile << "Netlist file: " << buffer2 << "   " << "Architecture file: data/arch/k4-n1.xml" << endl;
        //     ofile << "Array size: " << dim-2 << " x " << dim-2 << " logic blocks" << endl << endl;
        //     ofile << "#block name\t" << "x\t" << "y\t" << "subblk\t" << "block number\n";
        //     ofile << "#----------\t" << "--\t" << "--\t" << "------\t" << "------------\n";
        //     for(int i=0; i<nodes; i++){
        //         ofile << names[i] << "\t\t" << positions[idxMinCost][i]%dim << "\t" << positions[idxMinCost][i]/dim << "\t" << "0" << "\t\t#" << i << "\n";
        //     }             
        // } else{
        //     // printGrid(0);
        //     ofile1.open(buffer1);
        //     ofile1 << "Netlist file: " << buffer2 << "   " << "Architecture file: data/arch/k4-n1.xml" << endl;
        //     ofile1 << "Array size: " << dim-2 << " x " << dim-2 << " logic blocks" << endl << endl;
        //     ofile1 << "#block name\t" << "x\t" << "y\t" << "subblk\t" << "block number\n";
        //     ofile1 << "#----------\t" << "--\t" << "--\t" << "------\t" << "------------\n";
        //     for(int i=0; i<nodes; i++){
        //         ofile1 << names[i] << "\t\t" << positions[idxMinCost][i]%dim << "\t" << positions[idxMinCost][i]/dim << "\t" << "0" << "\t\t#" << i << "\n";
        //     } 
        // }
        ofile.close();
        ofile1.close();
    }

    delete h_edgeA;
    delete h_edgeB;
    delete v;
    delete v_i;

    printf("%s,%d,%d,%d,%lf,%.2lf,%ld\n", argv[3], nodes, edges, gridSize, double(cost_min[1])/edges, time_total[1], numberSwaps[1]);
    return 0;
}