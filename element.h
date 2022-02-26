#include"node.h"


struct element
{
    int value, nn;
    node *nodes;
    double **H;
    double **Hbc;
    double **C;
    double *P;
    element(){}
    element(int value, int b, int c, node *arr, int nn);
    ~element();
    void show();
    void showMore();
    void printH();
    void printHbc();
    void addMatrixes();
    void aggregate(double **arr);
    void CalculateP(double *n,double *n2, double *n3, double t, int surfaceId,int pc);
    double getLength(int surfaceId);
    void aggregateVector(double *vec);
    void showAggregationVector();
    void aggregateC(double **arr);
    void calculateC(int i, double **n, double c, double ro, double **J,int pc, double w);
};

element::element(int value, int b, int c, node *arr, int nn)
{
    
    this->value = value;
    this->nn = nn;
    nodes = new node[4];
    nodes[0] = arr[b-1];
    nodes[1] = arr[b+c-1];
    nodes[2] = arr[nodes[1].id];
    nodes[3] = arr[b];

    H = new double*[4];
    Hbc = new double*[4];
    C = new double*[4];
    P = new double[4];
    for(int i=0;i<4;++i)
    {
        H[i] = new double[4];
        Hbc[i] = new double[4];
        C[i] = new double[4];
        P[i] = 0;
        for(int j=0;j<4;++j)
        {
            H[i][j] = 0;
            Hbc[i][j] = 0;
            C[i][j] = 0;
        }
    }
}

element::~element()
{
    // delete[](nodes);
    // for(int i=0;i<4;++i)
    // {
    //     delete[](H[i]);
    //     delete[](Hbc[i]);
    // }
    // delete[]H;
    // delete[]Hbc;
}

void element::show()
{
    cout<<"value: "<<value<<endl;
    cout<<"nodes: "<<nodes[0].id<<", "<<nodes[1].id<<", "<<nodes[2].id<<", "<<nodes[3].id<<endl;    
}

void element::showMore()
{
    element::show();
    nodes[0].show();
    cout<<"\n";
    nodes[1].show();
    cout<<"\n";
    nodes[2].show();
    cout<<"\n";
    nodes[3].show();
    cout<<"\n";
}


void element::printH()
{
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
            cout<<H[i][j]<<"\t";
        cout<<endl;
    }
    cout<<endl;
}

void element::printHbc()
{
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
            cout<<Hbc[i][j]<<" ";
        cout<<endl;
    }
    cout<<endl;
}


void element::addMatrixes()
{
    for(int i=0;i<4;++i)
        for(int j=0;j<4;++j)
            Hbc[i][j] += H[i][j];
}

void element::aggregate(double **arr)
{
    addMatrixes();
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
        {
            int x = nodes[i].id-1;
            int y = nodes[j].id-1;
            arr[y][x] = Hbc[i][j];
        }
    }
}



double element::getLength(int surfaceId) 
{
    int firstNodeIndex=0, secondNodeIndex=0;

    switch (surfaceId) 
    {
        //bottom
        case 0: 
        {
            firstNodeIndex = 0;
            secondNodeIndex = 1;
            break;
        }
        //right
        case 1: 
        {
            firstNodeIndex = 1;
            secondNodeIndex = 2;
            break;
        }
        //up
        case 2: 
        {
            firstNodeIndex = 2;
            secondNodeIndex = 3;
            break;
        }
        //left
        case 3: 
        {
            firstNodeIndex = 3;
            secondNodeIndex = 0;
            break;
        }
    }

    double xDiff = nodes[firstNodeIndex].x - nodes[secondNodeIndex].x;
    double yDiff = nodes[firstNodeIndex].y - nodes[secondNodeIndex].y;
    return sqrt(xDiff * xDiff + yDiff * yDiff);
}

void element::CalculateP(double *n1, double *n2, double *n3, double t, int surfaceId,int pc)
{
    double det = getLength(surfaceId)/2;
    if(pc == 2)
    {
        double w1 = 1, w2 = 1;
        
        for(int i=0;i<4;++i)
        {
            P[i] += 300.0 * ( (w1 * n1[i]) * t + (w2 * n2[i]) * t ) * det;
        }
    }
    else if(pc == 3)
    {
        double w1 = 5.0/9.0, w2 = 8.0/9.0, w3 = w1;
        for(int i=0;i<4;++i)
        {
            P[i] += 300.0 * ( (w1 * n1[i]) * t + (w2 * n2[i]) * t + (w1 * n3[i]) * t) * det;

        }
    }
}

void element::aggregateVector(double *vec)
{
    for(int i=0;i<4;++i)
        vec[nodes[i].id-1] = P[i];
}

void element::showAggregationVector()
{
    for(int i=0;i<4;++i)
        cout<<P[i]<<" ";
    cout<<endl<<endl;
}

void element::aggregateC(double **arr)
{
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
        {
            int x = nodes[i].id-1;
            int y = nodes[j].id-1;
            arr[y][x] = C[i][j];
        }
    }
}

void element::calculateC(int i,double **n, double c, double ro, double **J, int pc, double w)
{
    double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if(pc == 2)
    {
        for(int j=0;j<4;++j)
        {
            for(int k=0;k<4;++k)
            {
                C[j][k] += c*ro*det*(n[i][j] * n[i][k]);
            }
        }
    }
    else if(pc == 3)
    {
        for(int j=0;j<4;++j)
        {
            for(int k=0;k<4;++k)
            {
                C[j][k] += c*ro*det*(n[i][j] * n[i][k]) * w;
            }
        }
    }

}