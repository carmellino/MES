#include"element.h"
#include <iomanip>  

struct grid
{   
    double initialTemp, simulationStepTime, ambientTemp, alfa, specificHeat, conductivity, density;
    double h, b;
    int nh, nb, nn, ne, pc;
    node *nodes;
    element *elements;
    double **nc;
    double **aggregationArr;
    double *aggregationVector;
    double *finalP;
    double **aggregationArrC;
    double *T0;

    grid(double H, double B, int NH, int NB, int pc);
    grid(double H, double B, int NH, int NB,int pc, double initialTemp, double simulationStepTime, double ambientTemp, double alfa, double specificHeat, double conductivity, double density);
    ~grid();
    void show();
    void showBcGrid();
    void calculateJakobian(int i, int j, double **dEta, double** dKsi,double **J);
    void calculateJakobians(double **dEta, double** dKsi);
    void inverseJakobian(double **J, double **JIn);
    void calculateDxDy(double **dx, double **dy, double **poEta, double** poKsi, double **JIn, int i);
    void calculateH(double **h, double **dx, double **dy, double **J,int a,double w);
    void multiplyH(double **h, double **J, double t);
    void aggregateAll();
    void addToAggregationArr(double **arr);
    void printAggregationArr();
    void printAggregationVector();
    void addToAggregationVector(double *vec);
    void addToAggregationArrC(double **arr);
    void printAggregationArrC();
    void addCToH();
    void pPlusCTimesT0();
    void aggregateAllC();
    void printFinalP();
    
};

grid::grid(double H, double B, int NH, int NB,int pc)
{
    this->pc = pc;
    h = H;
    b = B;
    nh = NH;
    nb = NB;
    nn = nh * nb;
    ne = (nh-1) * (nb-1);
    nodes = new node[nn];
    elements = new element[ne];
    
    //filling nodes 
    double dx = b/(nb-1), dy = h/(nh-1);
    int k=0;
    for(int i = 0; i < nb; ++i)
        for(int j=0; j<nh; ++j)
        {
            nodes[k] = node(k+1, dx * i, dy * j, h, b);
            k++;
        }

    //filling elements
    int ld = 1;
    for (int i = 0; i < ne; ++i)
    {
        if(i!= 0)
            if(i%(nb-1)  == 0)
                ld++;
        elements[i] = element(i+1, ld, nh, nodes,nn);
        ld = elements[i].nodes[3].id;
    }

    aggregationVector = new double[nn];
    finalP = new double[nn];
    aggregationArr = new double*[nn];
    aggregationArrC = new double*[nn];
    T0 = new double[nn];
    for(int i=0;i<nn;++i)
    {   
        aggregationVector[i] = 0.0;
        finalP[i] = 0.0;
        T0[i] = 100.0;
        aggregationArr[i] = new double[nn];
        aggregationArrC[i] = new double[nn];
        for(int j=0;j<nn;++j)
        {
            aggregationArr[i][j] = 0.0;
            aggregationArrC[i][j] = 0.0;
        }
    }
}

grid::grid(double H, double B, int NH, int NB,int pc, double initialTemp, double simulationStepTime, double ambientTemp, double alfa, double specificHeat, double conductivity, double density)
{
    this->pc = pc;
    this->initialTemp = initialTemp;
    this->simulationStepTime = simulationStepTime;
    this->ambientTemp = ambientTemp;
    this->alfa = alfa;
    this->specificHeat = specificHeat;
    this->conductivity = conductivity;
    this->density = density;
    
    h = H;
    b = B;
    nh = NH;
    nb = NB;
    nn = nh * nb;
    ne = (nh-1) * (nb-1);
    nodes = new node[nn];
    elements = new element[ne];
    
    //filling nodes 
    double dx = b/(nb-1), dy = h/(nh-1);
    int k=0;
    for(int i = 0; i < nb; ++i)
        for(int j=0; j<nh; ++j)
        {
            nodes[k] = node(k+1, dx * i, dy * j, h, b);
            k++;
        }

    //filling elements
    int ld = 1;
    for (int i = 0; i < ne; ++i)
    {
        if(i!= 0)
            if(i%(nb-1)  == 0)
                ld++;
        elements[i] = element(i+1, ld, nh, nodes,nn);
        ld = elements[i].nodes[3].id;
    }

    aggregationVector = new double[nn];
    finalP = new double[nn];
    aggregationArr = new double*[nn];
    aggregationArrC = new double*[nn];
    T0 = new double[nn];
    for(int i=0;i<nn;++i)
    {   
        aggregationVector[i] = 0.0;
        finalP[i] = 0.0;
        T0[i] = initialTemp;
        aggregationArr[i] = new double[nn];
        aggregationArrC[i] = new double[nn];
        for(int j=0;j<nn;++j)
        {
            aggregationArr[i][j] = 0.0;
            aggregationArrC[i][j] = 0.0;
        }
    }
}

grid::~grid()
{
    //delete[]nodes;
    //delete[]elements;
}

void grid::show()
{
    for(int i=0; i<ne; ++i)
    {
        elements[i].showMore();
        cout<<endl;
    }
}

void grid::showBcGrid()
{
    cout<<"BC grid:\n";
    for(int i=nh; i>0; --i)
    {
        int k=i-1;
        for(int j=0;j<nb;++j)
        {
            cout << nodes[k].bc;
            k+=nh;
        }
        cout<<endl;
    }
    cout<<endl;
}

 void grid::calculateJakobian(int i, int j, double **dEta, double** dKsi,double **J)
{
    J[0][0] = 0;
    J[1][0] = 0;
    J[0][1] = 0;
    J[1][1] = 0;
 
    for(int k=0;k<4;++k)
    {
        J[0][0] += dKsi[j][k] * elements[i].nodes[k].y;
        J[0][1] += dEta[j][k] * elements[i].nodes[k].y;
        J[1][0] += dKsi[j][k] * elements[i].nodes[k].x;
        J[1][1] += dEta[j][k] * elements[i].nodes[k].x;
    }
}

void grid::inverseJakobian(double **J, double **JIn)
{
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    det = 1/det;
    JIn[0][0] = det * J[0][0];
    JIn[0][1] = det * J[0][1];
    JIn[1][0] = det * J[1][0];
    JIn[1][1] = det * J[1][1];
}


void grid::calculateJakobians(double **dEta, double** dKsi)
{
    double **J = new double*[2];
    double **JIn = new double*[2];
    double **dx = new double*[pc*pc];
    double **dy = new double*[pc*pc];
    J[0] = new double[2];
    J[1] = new double[2];
    JIn[0] = new double[2];
    JIn[1] = new double[2];
    for(int i =0;i<pc*pc; ++i)
    {
        dx[i] = new double[4];
        dy[i] = new double[4];
    }
    
    for(int i=0;i<ne;++i)
    {
        double w1 = 5.0/9.0, w2 = 8.0/9.0, w3 = w1, w4 = w1;
        for(int j=0;j<pc*pc; ++j)
        {
            calculateJakobian(i,j,dEta,dKsi,J);
            inverseJakobian(J,JIn);
            // cout<<"("<<i<<", "<<j<<")"<<endl;
            //cout<<J[0][0] << "\t"<<J[0][1]<<endl<<J[1][0]<<"\t"<<J[1][1]<<endl<<endl;
            //cout<<JIn[0][0] << "\t"<<JIn[0][1]<<endl<<JIn[1][0]<<"\t"<<JIn[1][1]<<endl<<endl<<endl;
            calculateDxDy(dx,dy,dEta,dKsi,JIn,j);
            calculateH(elements[i].H,dx,dy,J, j,w1 * w4);
            
            elements[i].calculateC(j,nc,specificHeat,density,J,pc, w1 * w4);
            w1 = w2;
            w2 = w3;
            w3 = w1;
            if(j == 7)
            {
                w1 = 8.0/9.0;
                w4 = w1;
            }
        }

        multiplyH(elements[i].H,J,conductivity);
        
        elements[i].printH();
    }

    delete[]J[0];
    delete[]J[1];
    delete[]JIn[0];
    delete[]JIn[1];
    for(int i=0;i<pc*pc;++i)
    {
        delete[]dx[i];
        delete[]dy[i];
    }
    delete[]dx;
    delete[]dy;
    delete[]J;
    delete[]JIn;
}

void grid::calculateDxDy(double **dx, double **dy, double **dEta, double** dKsi, double **JIn, int i)
{
    for(int j=0;j<4;++j)
    {
        dx[i][j] = JIn[0][0] * dEta[i][j] + JIn[0][1] * dKsi[i][j];
        dy[i][j] = JIn[1][0] * dEta[i][j] + JIn[1][1] * dKsi[i][j];
    }
}

void grid::calculateH(double **h, double **dx, double **dy, double **J, int a,double w)
{   
    if(pc == 2)
    {
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<4; ++j)
            {
                h[i][j] += (dx[a][i] * dx[a][j]);
                h[i][j] += (dy[a][i] * dy[a][j]);
            }
        }
    }
    else if(pc == 3)
    {
        for(int i=0; i<4; ++i)
        {
            for(int j=0; j<4; ++j)
            {
                h[i][j] += (dx[a][i] * dx[a][j]) * w;
                h[i][j] += (dy[a][i] * dy[a][j]) * w;
            }
        }
    }
}

void grid::multiplyH(double **h, double **J, double t)
{
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    for(int i=0; i<4; ++i)
    {
        for(int j=0; j<4; ++j)
        {
            h[i][j]*=det;
            h[i][j]*=t;
        }
    }
}

void grid::aggregateAll()
{
    for(int i=0; i<ne; ++i)
    {
        double **aggrArr = new double*[nn];
        for(int i=0;i<nn;++i)
        {
            aggrArr[i] = new double[nn];
            for(int j=0;j<nn;++j)
                aggrArr[i][j] = 0;
        }
        elements[i].aggregate(aggrArr);
        addToAggregationArr(aggrArr);

        for(int i=0;i<nn;++i)
            delete[]aggrArr[i];
        delete[]aggrArr;
    }
    //printAggregationArr();
}

void grid::addToAggregationArr(double **arr)
{
    for(int j=0;j<nn;++j)
        for(int k=0;k<nn;++k)
            aggregationArr[j][k] += arr[j][k];
}

void grid::printAggregationArr()
{
    for(int j=0;j<nn;++j)
    {
        for(int k=0;k<nn;++k)
            cout<<setw(9)<<aggregationArr[j][k]<<" ";
        cout<<endl;
    }
    cout<<endl;
}


void grid::addToAggregationVector(double *vec)
{
    for(int j=0;j<nn;++j)
        aggregationVector[j] += vec[j];
}

void grid::printAggregationVector()
{
    for(int j=0;j<nn;++j)
        cout<<aggregationVector[j]<<" ";
    cout<<endl<<endl;
}

void grid::aggregateAllC()
{
    for(int i=0; i<ne; ++i)
    {
        double **aggrArr = new double*[nn];
        for(int i=0;i<nn;++i)
        {
            aggrArr[i] = new double[nn];
            for(int j=0;j<nn;++j)
                aggrArr[i][j] = 0;
        }

        elements[i].aggregateC(aggrArr);
        addToAggregationArrC(aggrArr);

        for(int i=0;i<nn;++i)
            delete[]aggrArr[i];
        delete[]aggrArr;
    }
    //printAggregationArrC();
}

void grid::addToAggregationArrC(double **arr)
{
    for(int j=0;j<nn;++j)
        for(int k=0;k<nn;++k)
            aggregationArrC[j][k] += arr[j][k];
}

void grid::printAggregationArrC()
{
    for(int j=0;j<nn;++j)
    {
        for(int k=0;k<nn;++k)
            cout<<setw(8)<<aggregationArrC[j][k]<<" ";
        cout<<endl;
    }
    cout<<endl;
}

void grid::addCToH()
{
    for(int j=0;j<nn;++j)
        for(int k=0;k<nn;++k)
            aggregationArr[j][k] += aggregationArrC[j][k]/simulationStepTime;
    //printAggregationArr();
}

void grid::pPlusCTimesT0()
{
    double sum;

    for(int i=0;i<nn;++i)
    {
        sum = 0.0;
        for(int j=0;j<nn;++j)
        {
            sum += (aggregationArrC[i][j]/simulationStepTime) * T0[j];
        }
        finalP[i] = aggregationVector[i] + sum;
    }
    
}

void grid::printFinalP()
{
    for(int j=0;j<nn;++j)
        cout<<finalP[j]<<" ";
    cout<<endl<<endl;
}