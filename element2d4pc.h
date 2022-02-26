#include<iostream>
#include"math.h"
#include"grid.h"
using namespace std;

struct wall
{
    element *w;
    double *n1;
    double *n2;
    double *n3;
    int size;

    wall(){};
    wall(int size, double ksi1, double eta1, double ksi2, double eta2);
    wall(int size, double ksi1, double eta1, double ksi2, double eta2, double ksi3, double eta3);

};

struct element2d4pc
{
    double **dEta;
    double **dKsi; 
    double **nc;
    
    wall lWall;
    wall bWall;
    wall rWall;
    wall tWall;
    grid *g;


    element2d4pc(grid *g);
    ~element2d4pc();
    double fun(double a);
    double fun1(double a);
    double fun3(double a);
    void calculateDEtaAndDKsi(int pc);
    void printDEtaAndDKsi();
    void fillWalls(element *e);
    void calculateHbc(element e,wall w,int surfaceId);
    void calculateAllHbc();
    void aggregateAllVectors();
    void calculateAllVectorsP();
    void fillN(double **n, double ksi, double eta);
    void s();
    double* equationSolve(double **a, double *b,int num);
    
    
};

wall::wall(int size, double ksi1, double eta1, double ksi2, double eta2)
{
    this->size = size;
    w = new element[size];
    n1 = new double[4];
    n2 = new double[4];
    
    n1[0] = 0.25 * (1.0-ksi1) * (1-eta1);
    n1[1] = 0.25 * (1+ksi1) * (1-eta1);
    n1[2] = 0.25 * (1+ksi1) * (1+eta1);
    n1[3] = 0.25 * (1-ksi1) * (1+eta1);
    
    n2[0] = 0.25 * (1.0-ksi2) * (1-eta2);
    n2[1] = 0.25 * (1+ksi2) * (1-eta2);
    n2[2] = 0.25 * (1+ksi2) * (1+eta2);
    n2[3] = 0.25 * (1-ksi2) * (1+eta2);
    // cout<<"\nN\n";
    // cout<<"ksi "<<ksi1<<" eta "<<eta1<<endl;
    // cout<<n1[0]<<" "<<n1[1]<<" "<<n1[2]<<" "<<n1[3]<<endl<<endl;
    // cout<<"ksi2 "<<ksi2<<" eta2 "<<eta2<<endl;
    // cout<<n2[0]<<" "<<n2[1]<<" "<<n2[2]<<" "<<n2[3]<<endl<<endl;

}

wall::wall(int size, double ksi1, double eta1, double ksi2, double eta2, double ksi3, double eta3)
{
    this->size = size;
    w = new element[size];
    n1 = new double[4];
    n2 = new double[4];
    n3 = new double[4];
    
    n1[0] = 0.25 * (1.0-ksi1) * (1-eta1);
    n1[1] = 0.25 * (1+ksi1) * (1-eta1);
    n1[2] = 0.25 * (1+ksi1) * (1+eta1);
    n1[3] = 0.25 * (1-ksi1) * (1+eta1);
    
    n2[0] = 0.25 * (1.0-ksi2) * (1-eta2);
    n2[1] = 0.25 * (1+ksi2) * (1-eta2);
    n2[2] = 0.25 * (1+ksi2) * (1+eta2);
    n2[3] = 0.25 * (1-ksi2) * (1+eta2);

    n3[0] = 0.25 * (1.0-ksi3) * (1-eta3);
    n3[1] = 0.25 * (1+ksi3) * (1-eta3);
    n3[2] = 0.25 * (1+ksi3) * (1+eta3);
    n3[3] = 0.25 * (1-ksi3) * (1+eta3);

    // cout<<"\nN\n";
    // cout<<"ksi "<<ksi1<<" eta "<<eta1<<endl;
    // cout<<n1[0]<<" "<<n1[1]<<" "<<n1[2]<<" "<<n1[3]<<endl<<endl;
    // cout<<"ksi2 "<<ksi2<<" eta2 "<<eta2<<endl;
    // cout<<n2[0]<<" "<<n2[1]<<" "<<n2[2]<<" "<<n2[3]<<endl<<endl;

}

element2d4pc::element2d4pc(grid *g)
{
    this->g = g;
    dKsi = new double*[g->pc*g->pc];
    dEta = new double*[g->pc*g->pc];
    nc = new double*[g->pc*g->pc];

    for(int i =0;i<g->pc*g->pc; ++i)
    {
        dKsi[i] = new double[4];
        dEta[i] = new double[4];
        nc[i] = new double[4];
    }

    if(g->pc == 2)
        fillN(nc,-(1.0/sqrt(3.0)),-(1.0/sqrt(3.0)));
    else if(g->pc == 3)
        fillN(nc,sqrt(3.0/5.0), 1.0);
    g->nc = nc;

    if(g->pc == 2)
    {
        bWall =  wall(g->nb, -(1.0/sqrt(3.0)), -1.0, (1.0/sqrt(3.0)),-1.0);
        rWall =  wall(g->nh, 1.0, -(1.0/sqrt(3.0)), 1.0,(1.0/sqrt(3.0)));
        tWall =  wall(g->nb, -(1.0/sqrt(3.0)), 1.0, (1.0/sqrt(3.0)),1.0);
        lWall =  wall(g->nh, -1.0, (1.0/sqrt(3.0)), -1.0, -(1.0/sqrt(3.0)));
    }
    else if(g->pc == 3)
    {
        bWall =  wall(g->nb, -(sqrt(3.0/5.0)), -1.0, 0, -1.0, (sqrt(3.0/5.0)), -1.0);
        rWall =  wall(g->nh, 1.0, -(sqrt(3.0/5.0)), 1.0, 0, 1.0,(sqrt(3.0/5.0)));
        tWall =  wall(g->nb, (sqrt(3.0/5.0)), 1.0, 0, 1.0, -(sqrt(3.0/5.0)), 1.0);
        lWall =  wall(g->nh, -1.0, (sqrt(3.0/5.0)), -1.0, 0, -1.0, -(sqrt(3.0/5.0)));
    }

    fillWalls(g->elements);
}

element2d4pc::~element2d4pc()
{
    for(int i =0;i<4; ++i)
    {
        delete[]dKsi[i];
        delete[]dEta[i];
        delete[]nc[i];
    }
    delete[]dKsi;
    delete[]dEta;
    delete[]nc;
}
double element2d4pc::fun(double a)
{
    return (1.0/4.0) * (1-a);
}

double element2d4pc::fun1(double a)
{
    return (1.0/4.0) * (1+a);
}

double element2d4pc::fun3(double a)
{
    return (1.0/4.0) * (a-1);
}

void element2d4pc::calculateDEtaAndDKsi(int pc)
{   
    
    if(pc == 2)
    {
        double x = 1.0/sqrt(3.0);
        double xyArr[4][2] = {{-x,-x},{x,-x}, {x,x}, {-x, x}};

        for (int i = 0; i < 4; i++)
        {
            dKsi[i][0] = fun3(xyArr[i][0]);
            dKsi[i][1] = -fun1(xyArr[i][0]);
            dKsi[i][2] = fun1(xyArr[i][0]);
            dKsi[i][3] = fun(xyArr[i][0]);

            dEta[i][0] = fun3(xyArr[i][1]);
            dEta[i][1] = fun(xyArr[i][1]);
            dEta[i][2] = fun1(xyArr[i][1]);
            dEta[i][3] = -fun1(xyArr[i][1]);
        }
    }
    else if(pc == 3)
    {
        double x = sqrt(3.0/5.0), x1 = 0.0;
        double xyArr[9][2] = {{-x,-x},{x1,-x}, {x,-x}, {x, x1}, {x,x}, {x1,x}, {-x,x}, {-x,x1}, {x1,x1}};
        for (int i = 0; i < 9; i++)
        {
            dKsi[i][0] = fun3(xyArr[i][0]);
            dKsi[i][1] = -fun1(xyArr[i][0]);
            dKsi[i][2] = fun1(xyArr[i][0]);
            dKsi[i][3] = fun(xyArr[i][0]);

            dEta[i][0] = fun3(xyArr[i][1]);
            dEta[i][1] = fun(xyArr[i][1]);
            dEta[i][2] = fun1(xyArr[i][1]);
            dEta[i][3] = -fun1(xyArr[i][1]);
        }
    }
}

void element2d4pc::printDEtaAndDKsi()
{
    if(g->pc == 2)
    {
        cout<<"\ndKsi:\n";
        for(int i=0; i<4; ++i)
        {
            for(int j=0;j<4;++j)
                cout<<dKsi[i][j]<<" ";
            cout<<endl;
        }
        cout<<"\ndEta:\n";
        for(int i=0; i<4; ++i)
        {
            for(int j=0;j<4;++j)
                cout<<dEta[i][j]<<" ";
            cout<<endl;
        }   
    }
    else if(g->pc == 3)
    {
        cout<<"\ndKsi:\n";
        for(int i=0; i<9; ++i)
        {
            for(int j=0;j<4;++j)
                cout<<dKsi[i][j]<<" ";
            cout<<endl;
        }
        cout<<"\ndEta:\n";
        for(int i=0; i<9; ++i)
        {
            for(int j=0;j<4;++j)
                cout<<dEta[i][j]<<" ";
            cout<<endl;
        }   
    }
}

void element2d4pc::fillWalls(element *e)
{
    int l = 0, r = 0, t=0,b=0;
    for(int i=0;i<g->ne;++i)
    {
        //bottom wall
        if((e[i].nodes[0].bc == 1) && (e[i].nodes[1].bc == 1))
        {
            bWall.w[b] = e[i];
            b++;
        }
        //right wall
        if((e[i].nodes[1].bc == 1) && (e[i].nodes[2].bc == 1))
        {
           rWall.w[r] = e[i];
            r++;
        }
        //top wall
        if((e[i].nodes[2].bc == 1) && (e[i].nodes[3].bc == 1))
        {
            tWall.w[t] = e[i];
            t++;
        }
        //left wall
        if((e[i].nodes[3].bc == 1) && (e[i].nodes[0].bc == 1))
        {
            lWall.w[l] = e[i];
            l++;
        }
    }
    // for(int i=0; i<g->nb-1; ++i)
    //     cout<<tWall.w[i].value<<" ";
    // cout<<endl;
    // for(int i=0; i<g->nb-1; ++i)
    //     cout<<bWall.w[i].value<<" ";
    // cout<<endl;
    // for(int i=0; i<g->nh-1; ++i)
    //     cout<<lWall.w[i].value<<" ";
    // cout<<endl;
    // for(int i=0; i<g->nh-1; ++i)
    //     cout<<rWall.w[i].value<<" ";
    // cout<<endl;
    
}


void element2d4pc::calculateHbc(element e, wall w, int surfaceId)
{
    double det = e.getLength(surfaceId);
    
    if(g->pc == 2)
    {
        for(int j=0;j<4;++j)
            for(int k=0;k<4;++k)
            {
                e.Hbc[j][k] += g->alfa * det/2 * w.n1[j] * w.n1[k];
                e.Hbc[j][k] += g->alfa * det/2 * w.n2[j] * w.n2[k];
            }
    }
    else if(g->pc == 3)
    {
        double w1 = 5.0/9.0, w2 = 8.0/9.0;
        for(int j=0;j<4;++j)
            for(int k=0;k<4;++k)
            {
                e.Hbc[j][k] += g->alfa * det/2 * w1 * w.n1[j] * w.n1[k];
                e.Hbc[j][k] += g->alfa * det/2 * w2 * w.n2[j] * w.n2[k];
                e.Hbc[j][k] += g->alfa * det/2 * w1 * w.n3[j] * w.n3[k];
            }
    }
}

void element2d4pc::calculateAllHbc()
{
    //bottom and top wall
    for(int i=0; i<g->nb-1; ++i)
    {
        calculateHbc(bWall.w[i],bWall,0);
        calculateHbc(tWall.w[i],tWall,2);
    }

    //left and right wall
    for(int i=0; i<g->nh-1; ++i)
    {
        calculateHbc(lWall.w[i],lWall,3);
        calculateHbc(rWall.w[i],rWall,1);
    }
    // for(int i=0;i<g->ne;++i)
    // {
    //     cout<<"HBC"<<endl;
    //     g->elements[i].printHbc();
    // }
    // //Wyświetlanie HBC
    // for(int i=0; i<g->nb-1; ++i)
    // {
    //     cout<<"HBC"<<endl;   
    //     bWall.w[i].printHbc();
    //     cout<<endl;
    //     cout<<"HBC"<<endl;   
    //     tWall.w[i].printHbc();
    //     cout<<endl;
    // }

    // //left and right wall
    // for(int i=0; i<g->nh-1; ++i)
    // {
    //     cout<<"HBC"<<endl;   
    //     lWall.w[i].printHbc();
    //     cout<<endl;
    //     cout<<"HBC"<<endl;   
    //     rWall.w[i].printHbc();
    //     cout<<endl;
    // }
    
}


void element2d4pc::aggregateAllVectors()
{
    calculateAllVectorsP();

    for(int i=0; i<g->ne; ++i)
    {
        double *aggrVec = new double[g->nn];
        for(int i=0;i<g->nn;++i)
           aggrVec[i] = 0;

        g->elements[i].aggregateVector(aggrVec);
        g->addToAggregationVector(aggrVec);
        delete[] aggrVec;
    }
    //g->printAggregationVector();
}

void element2d4pc::calculateAllVectorsP()
{
    //bottom and top wall
    for(int i=0; i<g->nb-1; ++i)
    {
        bWall.w[i].CalculateP( bWall.n1,  bWall.n2, bWall.n3,g->ambientTemp,0, g->pc);
        tWall.w[i].CalculateP( tWall.n1,  tWall.n2, tWall.n3,g->ambientTemp,2, g->pc);
    }
    //left and right wall
    for(int i=0; i<g->nh-1; ++i)
    {
        lWall.w[i].CalculateP( lWall.n1,  lWall.n2, lWall.n3, g->ambientTemp,3, g->pc);
        rWall.w[i].CalculateP( rWall.n1,  rWall.n2, rWall.n3, g->ambientTemp,1, g->pc);
    }
}

void element2d4pc::fillN(double **n, double ksi, double eta)
{
    if(g->pc == 2)
    {
        n[0][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[0][1] = 0.25 * (1+ksi) * (1-eta);
        n[0][2] = 0.25 * (1+ksi) * (1+eta);
        n[0][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = -ksi;
        n[1][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[1][1] = 0.25 * (1+ksi) * (1-eta);
        n[1][2] = 0.25 * (1+ksi) * (1+eta);
        n[1][3] = 0.25 * (1-ksi) * (1+eta);

        eta = -eta;
        n[2][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[2][1] = 0.25 * (1+ksi) * (1-eta);
        n[2][2] = 0.25 * (1+ksi) * (1+eta);
        n[2][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = -ksi;
        n[3][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[3][1] = 0.25 * (1+ksi) * (1-eta);
        n[3][2] = 0.25 * (1+ksi) * (1+eta);
        n[3][3] = 0.25 * (1-ksi) * (1+eta);
    }
    else if(g->pc == 3)
    {
        double sq = ksi, one = eta;

        ksi = -sq;
        eta = -sq;
        n[0][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[0][1] = 0.25 * (1+ksi) * (1-eta);
        n[0][2] = 0.25 * (1+ksi) * (1+eta);
        n[0][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = 0;
        n[1][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[1][1] = 0.25 * (1+ksi) * (1-eta);
        n[1][2] = 0.25 * (1+ksi) * (1+eta);
        n[1][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = sq;
        n[2][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[2][1] = 0.25 * (1+ksi) * (1-eta);
        n[2][2] = 0.25 * (1+ksi) * (1+eta);
        n[2][3] = 0.25 * (1-ksi) * (1+eta);

        eta = 0;
        n[3][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[3][1] = 0.25 * (1+ksi) * (1-eta);
        n[3][2] = 0.25 * (1+ksi) * (1+eta);
        n[3][3] = 0.25 * (1-ksi) * (1+eta);

        eta = sq;
        n[4][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[4][1] = 0.25 * (1+ksi) * (1-eta);
        n[4][2] = 0.25 * (1+ksi) * (1+eta);
        n[4][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = 0;
        n[5][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[5][1] = 0.25 * (1+ksi) * (1-eta);
        n[5][2] = 0.25 * (1+ksi) * (1+eta);
        n[5][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = -sq;
        n[6][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[6][1] = 0.25 * (1+ksi) * (1-eta);
        n[6][2] = 0.25 * (1+ksi) * (1+eta);
        n[6][3] = 0.25 * (1-ksi) * (1+eta);

        eta = 0;
        n[7][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[7][1] = 0.25 * (1+ksi) * (1-eta);
        n[7][2] = 0.25 * (1+ksi) * (1+eta);
        n[7][3] = 0.25 * (1-ksi) * (1+eta);

        ksi = 0;
        n[8][0] = 0.25 * (1.0-ksi) * (1-eta);
        n[8][1] = 0.25 * (1+ksi) * (1-eta);
        n[8][2] = 0.25 * (1+ksi) * (1+eta);
        n[8][3] = 0.25 * (1-ksi) * (1+eta);
    }
}

void element2d4pc::s()
{

    double *sol = equationSolve(g->aggregationArr, g->finalP, g->nn);
    double min = sol[0], max = 0.0;
    for(int i=0;i<g->nn;++i)
    {
        g->T0[i] = sol[i];
        g->finalP[i] = 0;
        //cout<<sol[i]<<endl;
        if (sol[i] < min)
            min = sol[i];
        else if(sol[i] > max)
            max = sol[i];
    }
    cout<<"MIN: "<<min<<" MAX: "<<max<<endl;
    cout<<endl;
    
}


double* element2d4pc::equationSolve(double **a, double *b,int num)
{
	double *n = new double[num];
    double *x1= new double[num];
    double *x2 = new double[num];
	double **M = new double*[num];
    for(int i=0;i<num;++i)
    {
        M[i] = new double[num];
        for(int j=0;j<num;++j)
            M[i][j] = 0.0;
    }

	for (int i = 0; i < num; i++)
		n[i] = 1 / a[i][i];

	// Calculate M = -D^-1 (L + U)
	for (int i = 0; i < num; i++)
		for (int j = 0; j < num; j++)
			if (i == j)
				M[i][j] = 0.0;
			else
				M[i][j] = -(a[i][j] * n[i]);

	for (int k = 0; k < 100; k++) 
    {
		for (int i = 0; i < num; i++) 
        {
			x2[i] = n[i] * b[i];
			for (int j = 0; j < num; j++)
				x2[i] += M[i][j] * x1[j];
		}
		for (int i = 0; i < num; i++)
			x1[i] = x2[i];
	}
	return x1;
}