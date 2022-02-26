#include<iostream>
using namespace std;

struct node
{
    int id, bc = 0;
    double x, y;

    node(){}
    node(int id, double x, double y, double h, double b);
    void show();
};

node::node(int id, double x, double y, double h, double b)
{
    this->id = id;
    this->x = x;
    this->y = y;

    if(x == 0)
        bc = 1;
    else if(y == 0)
        bc = 1;
    else if(x == b)
        bc = 1;
    else if(y == h)
        bc = 1;  
}

void node::show()
{
    cout<<"id: "<<id<<endl;
    cout<<"bc: "<<bc<<endl;
    cout<<"( "<<x<<", "<<y<<" )\n";
}