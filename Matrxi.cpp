#include "Matrix.h"

class Matrix
{
	private:
		int n,m;
		double **M, **L, **U;
	public:
		Matrix(int nn, int mm)
		{
			int i;
			n = nn;
			m = mm;
			M = new double*[n];
			for (i = 0; i < n; i++)  { M[i] = new double[m]; }
		}
		void input();
		void show();
		void LU();
		void QR();
		void SVD();
		~Matrix()
		{
			int i;
			for (i = 0; i < n; i++) delete[] M[i];
		}
};


void Matrix::input()
{
	cout << "�������������:" << endl;
	cin >> n;
	cout<< "�������������:" << endl;
	cin >> m;
	int i = 0,j=0;
	cout << "�������������Ԫ�أ�" << endl;

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			cin >> M[i][j];
		}
}

void Matrix::show()
{
	int i = 0,j=0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			cout << M[i][j] << ' ';
		}
		cout << endl;
	}


}
