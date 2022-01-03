#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

const double e =2.7182818284;


struct Path{
    double weight;
    double prevDelta;
    double delta;
    int nid;    //neuron id in *n
};


struct Neuron{

double value;   //Local input
double bias;
double bd;
double pbd;
double gradient;
double out; //Output value
Path *p;    //path array 
int nc; //number of paths/connections belonging to *this
};




struct Layer{

    int nneurons;   //number of neurons in cluster/layer
    Neuron *n;   //Local Neuron array to reference
    Layer *neighbor;


    void Transfer(int nid)  //compute target Neuron's input and assign it to it

    {

        double valueOut=0;
        Neuron *temp;
        temp=&n[nid];

        //for each connection, send w*v to paired neuron
        for(int i=0; i<n[nid].nc; i++)
        {
            valueOut=temp->p[i].weight * temp->out;

            //neuron nid(as referenced by p[]), which is in the other layer, receives the value
            neighbor->n[temp->p[i].nid].value+=valueOut;

        }


    }

    void Initialize(int size)
    {

            nneurons=size;
            n=new Neuron[nneurons];

            for(int i=0; i<nneurons; i++)
            {
                n[i].value=0.0;
                n[i].bias=1.0;
                n[i].out=0.0;

            }

    }

    void FormConnections(Layer& nl)//with neighboring layer
    {
        neighbor=&nl;

        int nCon=neighbor->nneurons;

        for(int i=0; i<nneurons; i++)
        {

            n[i].nc=nCon;
            n[i].p=new Path[nCon];

            //neuron 'i' will link its paths to neurons in the other layer
            for(int ii=0; ii<n[i].nc; ii++)
            {
                n[i].p[ii].weight=1.0;
                n[i].p[ii].prevDelta=0.0;
                n[i].p[ii].nid=ii;

            }

        }

    }

};



class Brain{
  public:
double eta;
double alpha;
Layer   input,
        hidden,
        output;

double *target;

void GetInput(double* in){
    for(int i=0; i<input.nneurons; i++)
        input.n[i].value=in[i];

}

void GetDesiredOutput(double* t)
{
    target=t;
}

void Initialize(int inputsize, int hiddensize, int outputsize)
{
    input.Initialize(inputsize);
    hidden.Initialize(hiddensize);
    output.Initialize(outputsize);

    input.FormConnections(hidden);
    hidden.FormConnections(output);

}

void BP()
{

    //Calculate gradients for output
    for(int i=0; i<output.nneurons; i++)
    {output.n[i].gradient=(target[i] - output.n[i].out) * (1 - output.n[i].out) * (1 + output.n[i].out);}

    Neuron* temp;
    for(int i=0; i<hidden.nneurons; i++)
        {
            temp=&hidden.n[i];
            temp->gradient=0;
            //for each connection...
            for(int ii=0; ii<hidden.n[i].nc; ii++)
            {
                //temp's gradient gets values in the form w1*g2 + w2*g2 + ... + wn*gn,
                //where w is the weight of p that leads to output.n[i] from temp(hidden), and g
                //is the gradient of that output at p[CurrentConnection].nid

                temp->gradient+= temp->p[ii].weight * output.n[temp->p[ii].nid].gradient;
            }

            //Multiply the resultant sums with d/dx S(x)
            temp->gradient*= (temp->out)*(1-temp->out);


        }
    //---------------------------------------------------------------------------
        //Calculate delta

for(int i=0; i<input.nneurons; i++)
    {
        temp=&input.n[i];



        for(int ii=0; ii<input.n[i].nc; ii++)
        {
            temp->p[ii].delta=eta* hidden.n[temp->p[ii].nid].gradient* temp->out;
            temp->p[ii].weight+=temp->p[ii].prevDelta*alpha +temp->p[ii].delta;
            temp->p[ii].prevDelta=temp->p[ii].delta;



        }

    }

    for(int i=0; i<hidden.nneurons; i++)
    {
        temp=&hidden.n[i];
        temp->bd=eta*temp->gradient;

        temp->bias+=temp->bd+ alpha*temp->pbd;
        temp->pbd=temp->bd;

        for(int ii=0; ii<hidden.n[i].nc; ii++){
            temp->p[ii].delta=eta* output.n[temp->p[ii].nid].gradient* temp->out;
            temp->p[ii].weight+=temp->p[ii].prevDelta*alpha+ temp->p[ii].delta;
            temp->p[ii].prevDelta=temp->p[ii].delta;


        }

    }


        for(int i=0; i<output.nneurons; i++)
        {
            temp=&output.n[i];

            temp->bias=eta*temp->gradient;
        }
    Zero(hidden);
    Zero(output);
    }

    void Process()
    {
        for(int i=0; i<input.nneurons; i++)
        {   input.n[i].out=input.n[i].value;
            input.Transfer(i);//transfer each neuron data in input to hidden
        }
        for(int i=0; i<hidden.nneurons; i++)
        {

            hidden.n[i].out=Sigmoid(hidden.n[i].value + hidden.n[i].bias);

            hidden.Transfer(i);
        }

        for(int i=0; i<output.nneurons; i++)
        {

            output.n[i].out=HyperTan(output.n[i].value + output.n[i].bias);
            cout<<"Output "<<i<<": "<<output.n[i].out<<endl;
        }


    }

    void Zero(Layer &l){ for(int i=0; i<l.nneurons; i++) l.n[i].value=0.0;}
    void Randomize(Layer &l)
    {
        for(int i=0; i<l.nneurons; i++)
        {
            for(int ii=0; ii<l.n[i].nc; ii++)
            {
                l.n[i].p[ii].weight=rand()%100/10;

            }
        }
    }
    Brain(){eta=0.9; alpha=0.4;}

     double Sigmoid(double x)
  {
    if (x < -45.0) return 0.0;
    else if (x > 45.0) return 1.0;
    else return (1.0 / (1.0 + pow(e, -x)));
  }
  double HyperTan(double x)
  {
    if (x < -10.0) return -1.0;
    else if (x > 10.0) return 1.0;
    else return tanh(x);
  }
};


int main(void){
  Brain b;

  double data[]={1.0,2.0, 3.0};
  double ddata[]={-0.25,0.14};



  b.Initialize(3,4,2);

  b.GetDesiredOutput(ddata);
  b.GetInput(data);


  b.Process();
  b.BP();
}