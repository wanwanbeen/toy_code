/***********************************************************
 Function:
 Strassen_Solver for Matrix Multiplication
 
 Jie Yang 2014
 
 C = A*B
 a.txt: contains matrix A
 b.txt: contains matrix B
 ************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

// Define a new matrix for later on pointer operations
float ** MatNew(int n)
{
    float ** M;
    M=new float * [n]();
    for (int i=0;i<n;i++)
        M[i]=new float[n]();
    return M;
}

// Define Matrix Addition
float ** MatAdd(float ** A,float ** B,int n)

{
    float ** C=MatNew(n);
    int i,j;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            C[i][j]=A[i][j]+B[i][j];
    return C;
}

// Define Matrix Subtraction
float ** MatSub(float ** A,float ** B,int n)
{
    float ** C=MatNew(n);
    int i,j;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            C[i][j]=A[i][j]-B[i][j];
    return C;
    
}

// Define Strassen method for Matrix Multiply
float ** MatStrassen(float** A,float** B,int n)
{
    float ** C=MatNew(n);
    
    if (n==1)
        C[0][0]=A[0][0]*B[0][0];  // single element in the matrix can simply multiply
    
    else
    {
        int m=n/2; // dividen and conquer
        
        // define pointers  for elements in matrix A B, and C, these elements are also matrix.
        float **A00,**A01,**A10,**A11;
        float **B00,**B01,**B10,**B11;
        float **C00,**C01,**C10,**C11;
        // define Qs in strassen method
        float **Q1,**Q2,**Q3,**Q4,**Q5,**Q6,**Q7;
        
        A00=MatNew(m);
        A01=MatNew(m);
        A10=MatNew(m);
        A11=MatNew(m);
        
        B00=MatNew(m);
        B01=MatNew(m);
        B10=MatNew(m);
        B11=MatNew(m);
        
        C00=MatNew(m);
        C01=MatNew(m);
        C10=MatNew(m);
        C11=MatNew(m);
        
        Q1 =MatNew(m);
        Q2 =MatNew(m);
        Q3 =MatNew(m);
        Q4 =MatNew(m);
        Q5 =MatNew(m);
        Q6 =MatNew(m);
        Q7 =MatNew(m);
        
        for (int i=0;i<m;i++)
        {
            for (int j=0;j<m;j++)
            {
                A00[i][j]=A[i][j]; //top-left quarter of A
                A01[i][j]=A[i][j+m]; //top-right quarter of A
                A10[i][j]=A[i+m][j];  //lower-left quarter of A
                A11[i][j]=A[i+m][j+m];  //lower-right quarter of A
                
                B00[i][j]=B[i][j];  //top-left quarter of B
                B01[i][j]=B[i][j+m];  //top-right quarter of B
                B10[i][j]=B[i+m][j];  //lower-left quarter of B
                B11[i][j]=B[i+m][j+m];  //lower-right quarter of B
            }
        }
        
        Q1=MatStrassen(MatAdd(A00,A11,m),MatAdd(B00,B11,m),m); //Q1=(A00+A11)*(B00+B11)
        Q2=MatStrassen(MatAdd(A10,A11,m),B00,m); //Q2=(A10+A11)*B00
        Q3=MatStrassen(A00,MatSub(B01,B11,m),m); //Q3=A00*(B01-B11)
        Q4=MatStrassen(A11,MatSub(B10,B00,m),m); //Q4=A11*(B10-B00)
        Q5=MatStrassen(MatAdd(A00,A01,m),B11,m); //Q5=(A00+A01)*B11
        Q6=MatStrassen(MatSub(A10,A00,m),MatAdd(B00,B01,m),m); //Q6=(A10-A00)*(B00+B01)
        Q7=MatStrassen(MatSub(A01,A11,m),MatAdd(B10,B11,m),m); //Q7=(A01-A11)*(B10+B11)
        
        
        C00=MatAdd(MatSub(MatAdd(Q1,Q4,m),Q5,m),Q7,m); //C00=Q1+Q4-Q5+Q7
        C10=MatAdd(Q2,Q4,m); //C10=Q2+Q4
        C01=MatAdd(Q3,Q5,m); //C01=Q3+Q5
        C11=MatAdd(MatSub(MatAdd(Q1,Q3,m),Q2,m),Q6,m); //C11=Q1+Q3-Q2+Q6
        
        for (int i=0;i<m;i++)
        {
            for (int j=0;j<m;j++)
            {
                C[i][j]=C00[i][j];
                C[i][j+m]=C01[i][j];
                C[i+m][j]=C10[i][j];
                C[i+m][j+m]=C11[i][j];
            }
        }
        
    }
    
    return C;
    
}

int main()
{
    ifstream Afile;
    ifstream Bfile;
    
    int n=0;
    
    // count dimension n
    Afile.open("a.txt");
    string TempString;
    while (getline(Afile,TempString))
        n++;
    Afile.close();
    
    cout<<"Dim ="<<n<<endl;
    
    // read file a.txt
    Afile.open("a.txt");
    float ** A=MatNew(n);
    for (int i=0;i<n;i++)
    {
        string Aline;
        getline(Afile,Aline);
        stringstream Alinestream(Aline);
        
        for (int j=0;j<n;j++)
        {
            getline(Alinestream,TempString,',');
            A[i][j]=stof(TempString);
        }
    }
    Afile.close();
 
    //read file b.txt
     Bfile.open("b.txt");
    float ** B=MatNew(n);
    for (int i=0;i<n;i++)
    {
        string Bline;
        getline(Bfile,Bline);
        stringstream Blinestream(Bline);
        for (int j=0;j<n;j++)
        {
            getline(Blinestream,TempString,',');
            B[i][j]=stof(TempString);
        }
    }
    Bfile.close();
    
    //C=A*B using Strassen method
    float ** C=MatNew(n);
    C=MatStrassen(A,B,n);
    
    //display A in a.txt
    cout<<"Matrix A="<<setw(10)<<endl;
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            cout<<A[i][j]<<setw(10);
        }
        cout<<endl;
    }
    
    //display B in b.txt
    cout<<endl;
    cout<<"Matrix B="<<setw(10)<<endl;
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            cout<<B[i][j]<<setw(10);
        }
        cout<<endl;
        
    }
    
    //display C=A*B using Strassen method
    cout<<endl;
    cout<<"Matrix C="<<setw(10)<<endl;
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            cout<<C[i][j]<<setw(10);
        }
        cout<<endl;
        
    }
    return 0;
}
