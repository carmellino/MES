#include<iostream>

#include"element2d4pc.h"
using namespace std;



int main()
{
    //1 
    double initialTemp = 100, simulationTime = 500, simulationStepTime = 50, ambientTemp = 1200, alfa = 300, specificHeat = 700, conductivity = 25, density = 7800;
    double h = 0.1 ,b = 0.1;
    int nh = 4, nb = 4, pc = 2;

    //2
    // double initialTemp = 100.0, simulationTime = 20.0, simulationStepTime = 1.0, ambientTemp = 1200.0, alfa = 300.0, specificHeat = 700.0, conductivity = 25.0, density = 7800.0;
    // double h = 0.1 ,b = 0.1;
    // int nh = 31, nb = 31, pc = 2;

    // grid g(h,b,nh,nb);
    // // g.show();
    // // g.showBcGrid();
    // grid g1(0.025,0.025,2,2);
    
    cout<<"Grid constructor\n";
    grid g2(h, b, nh, nb, pc, initialTemp, simulationStepTime, ambientTemp, alfa, specificHeat, conductivity, density);
    //g2.show();

    cout<<"Universal element constructor\n";
    element2d4pc e(&g2);

    cout<<"calculating dEta and dKsi\n";
    e.calculateDEtaAndDKsi(pc);
    e.printDEtaAndDKsi();
        
    cout<<"calculating jakobians\n";
    g2.calculateJakobians(e.dEta,e.dKsi);

    cout<<"calculating Hbc\n";
    e.calculateAllHbc();

    cout<<"H aggregation\n";
    g2.aggregateAll();
    
    cout<<"P aggregation\n";
    e.aggregateAllVectors();

    cout<<"C aggregation\n";
    g2.aggregateAllC();
    g2.printAggregationArrC();

    cout<<"H + C/dT\n";
    g2.addCToH();
    g2.printAggregationArr();
    
    for(int i=0;i<simulationTime/simulationStepTime;++i)
    {   
        cout<<endl;
        cout<<"P + (C/dT * T0)\n";
        g2.pPlusCTimesT0();
        g2.printFinalP();

        cout<<"Temperatures\n";
        e.s();
    }
}








