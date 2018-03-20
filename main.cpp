#include"iostream"
#include<vector>
using namespace std;
#include "Matrix.h"


int main()
{
	int t = 4,n=0,m=0;
	while (1)
	{
		cout << "������ֽⷽ����1��LU��   2��Cholesky��   3��QR  4,Household�任 5��TwoDiag 6,On_Side_Jacopbi 0���˳�" << endl;
		cin >> t;
		if (!t)   break;
		cout << "����������ά��(n����������m��������)" << endl;
		cout << "n=";
		cin >> n;
		cout << "m=";
		cin >> m;
		Matrix* M = new Matrix(n,m);

		M->input();
		switch (t)
		{
		case 1:
			M->LU();
			break;
		case 2:
			M->Cholesky();
			break;
		case 3:
			M->QR();
			break;
		case 4:
		M->Householder().show();
		case 5:
		M->TwoDiag().show();
			break;
		case 6:
		Matrix T=M->On_Side_Jacobi()[1];
		}
		delete M;
	}

	system("pause");
}

