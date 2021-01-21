#include <vector>
#include <iostream>

using namespace std;

vector<vector<int>> tablemesh;
vector<vector<int>> table1hop;

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
            tablemesh[i][j] = distY+distX;
            table1hop[i][j] = distX/2+distX%2+distY/2+distY%2;
        }
    }
}

void print(vector<vector<int>> &A){
    for(int i=0; i<A.size(); i++){
        for(int j=0; j<A.size(); j++){
            if(j==A.size()-1) cout << A[i][j] << endl;
            else cout << A[i][j] << " ";
        }
    }
}

int main(){
    int dim = 10;
    buildMatrices(dim);
    print(table1hop);
    cout << endl;
    print(tablemesh);
    return 0;
}