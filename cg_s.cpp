#include <bits/stdc++.h>
#include <omp.h>
#include <cmath>


using namespace std;

const double NEARZERO = 1.0e-10; 
vector<vector<double> > matrix_mult(vector<vector<double> >a ,vector<vector<double> >b)
{
	int size1 = a.size();
	int size2 =b[0].size();
	int size3 = b.size();
	vector<vector<double> >d(size1,vector<double>(size2));
	int i,j,k;


	for(i=0;i<size1;i++)
		for( j=0;j<size2;j++)
			for(k=0;k<size3;k++)
				d[i][j]+=a[i][k]*b[k][j];
	
	return d;
}
vector<vector<double> > transpose(vector<vector<double> > v){
	int row = v.size();
	int col = v[0].size();
	vector<vector<double > >B(col,vector<double>(row));
	int i,j;

	for(i =0;i<row;i++){
		for( j=0;j<col;j++){
			B[j][i] = v[i][j];
		}
	}
	return B;
}
vector<vector<double > > matrix_add(vector<vector<double> >a,vector<vector<double> >b){
	int size = a.size();
	int size2 = a[0].size();
	int i =0,j=0;

	for(i=0;i<size;i++){
		for(j=0;j<size2;j++){
			a[i][j]+=b[i][j];
		}
	}
	return a;
}
vector<vector<double > > matrix_subs(vector<vector<double> >a,vector<vector<double> >b){
	int size = a.size();
	int size2 = a[0].size();
	int i =0,j=0;

	for(i=0;i<size;i++){
		for(j=0;j<size2;j++){
			a[i][j]-=b[i][j];
		}
	}
	return a;
}
vector<vector<double> > scalar_multi(int sca,vector<vector<double > > a){
	int size = a.size();
	int size2 = a[0].size();
	int i,j;

	for(i=0;i<size;i++){
		for(j=0;j<size2;j++){
			a[i][j] = sca*a[i][j];
		}
	}
	return a;
}
vector<vector<double> > scalar_div(vector<vector<double > > a,int sca){
	int size = a.size();
	int size2 = a[0].size();
	int i,j;

	for(i=0;i<size;i++){
		for(j=0;j<size2;j++){
			a[i][j] = a[i][j]/sca;
		}
	}
	return a;
}
double norm(vector<vector<double> >a){
	double sum=0.0;
	for(int i =0;i<a.size();i++){
		for(int j=0;j<a[0].size();j++){
			a[i][j] *=a[i][j];
			sum+=a[i][j];
		}
	}
	sum = sqrt(sum);
	return sum;
}
double innerproduct(vector<vector<double> >a,vector<vector<double> >b){
	double sum=0.0;
	for(int i =0;i<a.size();i++){
		for(int j=0;j<a[0].size();j++){
			a[i][j] *=b[i][j];
			sum+=a[i][j];
		}
	}
	return sum;
}
void c_g(vector< vector<double> > A, vector< vector<double> > B, double tol);

void print(vector<vector<double> > v){
	int size = v.size();
	int size2 = v[0].size();
	for(int i =0;i<size;i++){
		for(int j=0;j<size2;j++){
			printf("%f ",v[i][j]);
		}
		printf("\n");
	}
	return ;
}
int main()
{
	double n;	
	cin>>n;
double starttime;double endtime;
	vector< vector<double> > A(n,vector<double>(n));
	vector< vector<double> > B(n,vector<double>(1));
	for(double i =0;i<n;i++){
		for(double j=0;j<n;j++){
			cin>>A[i][j];
		}
	}
	for(double i =0;i<n;i++){
		cin>>B[i][0];
	}
	double tol =1.0e-10;
	starttime=omp_get_wtime();
	c_g(A,B,tol);
	endtime=omp_get_wtime();
	cout<<"Time Taken by serial  -  "<<endtime-starttime<<endl;
	return 0;
}



void c_g(vector< vector<double> > A, vector< vector<double> > B,double tol){
	int size = B.size();
	vector<vector<double> >X(size,vector<double>(1,0.0));
	vector<vector<double > > R = B;
	vector<vector<double > > P = R;
	vector<vector<double > >temp = matrix_mult(transpose(R),R);	
	double num = temp[0][0];
	for(int i =0;i<size;i++){
		vector<vector<double> >AP = matrix_mult(A,P);
		temp = matrix_mult(transpose(P),AP);
		double alpha = num/temp[0][0];
		X = matrix_add(X,scalar_multi(alpha,P));
		R = matrix_subs(R,scalar_multi(alpha,AP));
		temp = matrix_mult(transpose(R),R);
		if(sqrt(temp[0][0]) <tol){
			break;		
		}
		double beta = (temp[0][0]/num);
		P = matrix_add(R,scalar_multi(beta,P));
		num = temp[0][0];
	}
	print(X);
	return ;
}








