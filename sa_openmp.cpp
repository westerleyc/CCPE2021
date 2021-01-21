#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
#include <algorithm> 
#include <fstream>
#include <omp.h>

#define NGRIDS 1000

using namespace std;
using namespace std::chrono;

int results[NGRIDS][4];
double randomvec[1000000];

//Tabelas com os custos para mesh e hop
vector<vector<int>> tablemesh;
vector<vector<int>> table1hop;

int bestCost;
int idxMinCost;

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

void printGraphInfo(int nodes, int edges, int *h_edgeA, int *h_edgeB, vector<int> &A, int *v, int *v_i){
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
    printf("\n\n");
}

void printPlacementInfo(int nodes, int gridSize, int *grid, int *positions, int cost){
    printf("grid: ");
    for (int i = 0; i < gridSize; ++i) printf("%d ", grid[i]);
    printf("\n");
    printf("positions: ");
    for (int i = 0; i < nodes; ++i) printf("%d ", positions[i]);
    printf("\n");
    printf("cost: %d", cost);
    printf("\n\n");
}

int gridCost(int edges, int dim, int *h_edgeA, int *h_edgeB, int *positions, int arch){
    int cost_ = 0, increment=0, distManhattanI, distManhattanJ, chess;
    for(int k=0; k<edges; k++){
        int ifrom = positions[h_edgeA[k]]/dim; 
        int jfrom = positions[h_edgeA[k]]%dim; 
        int ito = positions[h_edgeB[k]]/dim; 
        int jto = positions[h_edgeB[k]]%dim;
        
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

void annealing2(int iresult, int nodes, int dim, int gridSize, int &cost, int *grid, int *positions, int *v_i, int *v, vector<int> &A, int arch){
    int *localGrid, *localPositions;
    int currentCost = cost, nextCost, swapCount=0, increment, distManhattanI, distManhattanJ, chess;
    localGrid = new int[gridSize];
    for(int i=0; i<gridSize; i++) localGrid[i] = grid[i];
    localPositions = new int[nodes];
    for(int i=0; i<nodes; i++) localPositions[i] = positions[i];    
    //random vector index
    double random, valor;
    int randomctrl = 0;
    double T=100;
    auto start = high_resolution_clock::now();
    while(T>=0.00001){
        for(int i=0; i<gridSize; i++){
            for(int j=i+1; j<gridSize; j++){                
                //if we're looking at 2 empty spaces, skip                   
                if(localGrid[i]==255 && localGrid[j]==255)
                    continue;

                int node1 = localGrid[i], node2 = localGrid[j];
                nextCost = currentCost;
                //remove cost from object edges                
                if(node1!=255){
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
                        nextCost -= increment; 
                    }
                }
                if(node2!=255){
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
                // cout << "debug6" << endl;
                //swap positions
                int old1, old2;
                old1 = i;
                old2 = j;
                if(node1!=255) localPositions[node1] = old2;
                if(node2!=255) localPositions[node2] = old1;
                localGrid[j] = node1;
                localGrid[i] = node2;
                //recalculate cost
                if(node1!=255){
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
                if(node2!=255){
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
                //parameter for annealing probability
                valor = exp(-1*(nextCost - currentCost)/T);
                //random number between 0 and 1
                random = randomvec[randomctrl];
                randomctrl++;
                if(randomctrl==1000000) randomctrl=0;

                //if cost after changes is less than before or if cost is higher but we're in the annealing probanility range, return
                if(nextCost <= currentCost || random <= valor){
                    currentCost = nextCost;
                    swapCount++;
                }
                //else, undo changes and stay with previous cost
                else{
                    if(node1!=255) localPositions[node1] = old1;
                    if(node2!=255) localPositions[node2] = old2;
                    localGrid[j] = node2;
                    localGrid[i] = node1;
                }
            }
            if(cost==0) break;             
            T*=0.999;
        }   
        //aqui T*=0.999;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    for(int i=0; i<gridSize; i++) grid[i] = localGrid[i];
    for(int i=0; i<nodes; i++) positions[i] = localPositions[i]; 
    
    results[iresult][1] = currentCost;
    results[iresult][2] = duration.count()/1000;
    results[iresult][3] = swapCount;

    if(currentCost<=bestCost){
        bestCost = currentCost;
        idxMinCost = iresult;
    }
    //if(iresult%100==0){
    // for(int i=0; i<gridSize; i++){
    //         cout << localGrid[i] << " ";
    //     } 
    // cout << endl;
    //}
    delete localPositions;
    delete localGrid;
}


int main(int argc, char** argv) {
    string arqName(argv[1]);
    srand (time(NULL));
    //graph variables 
    for(int i=0;i<1000000;i++){
        randomvec[i] = (double)rand() / (double)(RAND_MAX);
    }
    int nodes, edges;
    int *h_edgeA, *h_edgeB;
    vector<int> A;
    int *v, *v_i;

    //placement variables
    int dim, gridSize;
    //grids random simple
    int **grid;             
    int cost = 100000;   
    //READING INPUT
    if(argc<2){
        cout << "erro de entrada" << endl;
        return 0;
    }

    FILE *fptr;
	fptr = fopen(argv[1],"r");
    if (fptr == NULL) {
        printf("Error! opening file\n");
        exit(1);
    }
    //INPUT READ
    int c = 0, n1, n2;
    //NUMERO DE GRIDS
    int nGrids = NGRIDS;  
    while(fscanf(fptr, "%d %d", &n1, &n2) != EOF) {
        if (c == 0) {
		    nodes = n1;
            edges = n2;
            h_edgeA = new int[edges];
            h_edgeB = new int[edges];
            v = new int[nodes];
            v_i = new int[nodes];
            for(int i=0; i<nodes; i++){
                v[i]=0; v_i[i]=0;
            }
        } else {
            h_edgeA[c-1] = n1;
            h_edgeB[c-1] = n2;
            v[n1]++;
            if(n1!=n2) v[n2]++;
        }
      c++;
    }
    
    dim = ceil(sqrt(nodes))+1;
    buildMatrices(dim);
    gridSize = dim*dim;
    const int num_arch = 2;
    double time_total[num_arch] = {0.0, 0.0};
    int cost_min[num_arch] = {-1,-1};

    //Aloca memoria pros grids
    grid = new int*[nGrids];
    for(int i=0; i<nGrids; i++) grid[i] = new int[gridSize];

    bool *borders;   
    borders = new bool[gridSize];
    for(int i=0; i<gridSize; i++){
        borders[i] = false;
        if(i/dim==0) borders[i] = true;
        else if(i/dim==dim-1) borders[i] = true;
        if(i%dim==0) borders[i] = true;
        else if(i%dim==dim-1) borders[i] = true; 
    }
        

    for (int k = 0; k < num_arch; ++k) {
        //Preenche os grids
        for(int i=0; i<nGrids; i++){
            for(int j=0; j<gridSize; j++){
                if(j < nodes) grid[i][j] = j;
                else grid[i][j] = 255;
            }
            for (int j=gridSize-1; j>0; --j) swap (grid[i][j],grid[i][(j+rand()%gridSize)%gridSize]);
        }

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
        
        int positions[nGrids][nodes];
        for(int k=0; k<nGrids; k++){
            for(int i=0; i<gridSize; i++){
                for(int j=0; j<nodes; j++){
                    if(j==grid[k][i]) positions[k][j]=i;
                }
            }
        }

        //fclose(fptr);
        //END OF INPUT READ    

        ofstream ofile, ofile1;
        char buffer[strlen(argv[1])+20] = "\0";
        // char beg[10] = "results/"; 
        char out[7] = ".mesh";
        // strcat(buffer,beg);
        strcat(buffer,argv[1]);
        strcat(buffer,out);
        ofile.open(buffer);

        char buffer1[strlen(argv[1])+20] = "\0";
        // char beg1[10] = "results/"; 
        char out1[7] = ".1hop";
        // strcat(buffer,beg);
        strcat(buffer1,argv[1]);
        strcat(buffer1,out1);
        ofile1.open(buffer1);

        //ofile << "Number of Nodes: " << nodes << " Number of Edges: " << edges << endl;
        //for (int i = 0; i < edges; ++i) ofile << h_edgeA[i] << " -> " << h_edgeB[i] << endl;   
        //ofile << "v: ";
        //for (int i = 0; i < nodes; ++i) ofile << v[i] << " ";
        //ofile << endl;
        //ofile << "v_i: ";
        //for (int i = 0; i < nodes; ++i) ofile << v_i[i] << " ";
        //ofile << endl;
        //ofile << "A: ";
        //for (int i = 0; i < A.size(); ++i) ofile << A[i] << " ";
        //ofile << endl << endl;  

        bestCost=10000000;
        idxMinCost=0;
        auto start = high_resolution_clock::now();

        #pragma omp parallel for
        for (int i = 0; i < nGrids; ++i) {
            int cost = gridCost(edges, dim, h_edgeA, h_edgeB, positions[i], k);
            // printPlacementInfo(nodes, gridSize, grid[i], positions[i], cost);
            results[i][0] = cost;
            if(cost == 0){
                results[i][1] = cost;
                results[i][2] = 0;
                results[i][3] = 0;
                /*
                for(int j=0; j<gridSize; ++j){
                    cout << grid[i][j] << " ";
                } 
                cout << endl;
                */   
                continue;
            }
            annealing2(i, nodes, dim, gridSize, cost, grid[i], positions[i], v_i, v, A, k);
        }
        auto stop = high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = (stop-start);
        //auto duration = duration_cast<seconds>(stop - start);

        time_total[k] = duration.count()/1000;

        // ofile << "-----" << endl;
        // for(int i=0; i<nGrids; i++){
        //     ofile << gridCost(edges, dim, h_edgeA, h_edgeB, positions[i], k) << endl;
        // }
        // ofile << "-----" << endl;
        /*
        ofile << "Init\t" << "A2\t" << "t2\t" << "swap2\t" << endl;
        for(int i=0; i<nGrids; i++){
            for(int j=0; j<4; j++){
                ofile << results[i][j] << "\t";
            }
            ofile << endl;
        }
        */
        vector<int> costs;    
        int totalTime=0;
        for(int i=0; i<nGrids; i++){
            if(results[i][1]>=0) costs.push_back(results[i][1]);
            totalTime += results[i][2];
        }
        cost_min[k] = *min_element(costs.begin(), costs.end());
        // if(k==0){
        //     for(int i=0; i<nGrids; i++){
        //         for(int j=0; j<gridSize; j++){
        //             ofile << grid[i][j] << " ";
        //         }
        //         ofile << endl;
        //     } 
        // } else if(k==1){        
        //     for(int i=0; i<nGrids; i++){
        //         for(int j=0; j<gridSize; j++){
        //             ofile1 << grid[i][j] << " ";
        //         }
        //         ofile1 << endl;
        //     } 
        // }
        if(k==0){            
            ofile << "Netlist file: " << arqName << "   " << "Architecture File: " << endl;
            ofile << "Array size: " << dim-2 << " x " << dim-2 << " logic blocks" << endl << endl;
            ofile << "#block name\t" << "x\t" << "y\t" << "subblk\t" << "block number\n";
            ofile << "#----------\t" << "--\t" << "--\t" << "------\t" << "------------\n";
            for(int i=0; i<nodes; i++){
                ofile << i << "\t\t" << positions[idxMinCost][i]%dim << "\t" << positions[idxMinCost][i]/dim << "\n";
            }             
        } else{//} if(k==1){      
            ofile1 << "Netlist file: " << arqName << "   " << "Architecture File: " << endl;
            ofile1 << "Array size: " << dim-2 << " x " << dim-2 << " logic blocks" << endl << endl;
            ofile1 << "#block name\t" << "x\t" << "y\t" << "subblk\t" << "block number\n";
            ofile1 << "#----------\t" << "--\t" << "--\t" << "------\t" << "------------\n";
            for(int i=0; i<nodes; i++){
                ofile1 << i << "\t\t" << positions[idxMinCost][i]%dim << "\t" << positions[idxMinCost][i]/dim << "\n";
            } 
        }
        ofile.flush();
        ofile1.flush();

    }

    delete h_edgeA;
    delete h_edgeB;
    delete v;
    delete v_i;

    //ofile << "Menor custo: " << *min_element(costs.begin(), costs.end()) << endl;
    //ofile << "Tempo Total: " << totalTime << endl;
    //ofile.close();

    printf("%s,%d,%d,%.2lf,%.2lf\n", argv[2], cost_min[0], cost_min[1], time_total[0], time_total[1]);

    return 0;
}