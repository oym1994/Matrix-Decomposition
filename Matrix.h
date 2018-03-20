#pragma once
#include <iostream>
#include <fstream>
#include<cmath>
#include <time.h>
#include <windows.h>  

using namespace std;

int sign(double d)
{
	if (d<0) return -1;
	else if (d == 0) return 0;
	else return 1;
}

class Matrix
{
private:
	int n, m;  //n---行数，m------列数
	double **M, **L, **U,det;
	string name;
public:

	Matrix(int nn, int mm)
	{
		int i,j;
		n = nn;
		m = mm;
		M = new double*[n];
		for (i = 0; i < n; i++)  { M[i] = new double[m]; }
		for (i = 0; i < n; i++)
			for (j = 0; j < m;j++)
			{
				M[i][j] = 0;
			}
	}
	void input();
	void show();
	bool CouldLU();
	bool CouldQR();
	bool CouldCholesky();

	Matrix* LU();
	Matrix Trans();
	Matrix Givens(double x,double y);
	void ppp();
	void sss();
	void normalize();
	Matrix* On_Side_Jacobi();


	Matrix* QR();
	Matrix* FG();

	Matrix Cholesky();

    void SVD();
	void DET();
	Matrix Householder();      //Housedolder变换
	Matrix GetLine(int M);  //取矩阵第M列
	Matrix GetLine(int a,int b);  //取矩阵第a到第b列

	Matrix GetRow(int N ); //取矩阵第N行
	Matrix GetRow(int a, int b);  //取矩阵第a到第b列

	Matrix TwoDiag(); //矩阵二对角化


	double GetNorm();   //获得范数平方
    Matrix operator +(Matrix& M1);
	Matrix operator -(Matrix& M1);
	Matrix operator*(Matrix& M1);
	Matrix operator*(double d);

	Matrix operator/(double d);

	friend Matrix operator*(double z, Matrix& MM)//声明运算符函数为友元函数
	{
		Matrix TEMP(MM.n, MM.m);
		int i, j;
		for (i = 0; i < MM.n; i++)
			for (j = 0; j < MM.m; j++)
			{
				TEMP.M[i][j] = z*MM.M[i][j];
			}
		return TEMP;

	};

   /* Matrix operator*(Matrix& MM, double& z)
		{
		Matrix TEMP(MM.n, MM.m);
		int i, j;
		for (i = 0; i < MM.n; i++)
			for (j = 0; j < MM.m; j++)
			{
				TEMP.M[i][j] = z*MM.M[i][j];
			}
		return TEMP;
	     }
		 */

	void operator=(Matrix& M1);
	bool operator==(Matrix& M1);
	bool  IsPD();
	bool IsSymmetry();


	~Matrix()
	{
		int i;
		for (i = 0; i < n; i++) delete[] M[i];
	}
};

void Matrix::normalize()
{
	int i = 0, j = 0;
	double temp = 0;
	for (i = 0; i < m; i++)
	{
		temp = this->GetLine(i).GetNorm();
		for (j = 0; j < n; j++)
		{
			if (temp != 0)	M[j][i] = M[j][i] / temp ;
		}
	}
}


void Matrix::input()
{
	int i = 0, j = 0;
	cout << "下面请输入矩阵元素：" << endl;

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			cin >> M[i][j];
		}
}

void Matrix::show()
{
	int i = 0, j = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
		cout << M[i][j] << ' ';
		}
		cout << endl;
	}
}



Matrix Matrix::operator +(Matrix& M1)
{
	if (m != M1.m || n != M1.n)
	{
		cout << "两矩阵维数不一致，无法进行相加" << endl;
		exit(1);
	}
	int i = 0, j = 0;
	Matrix outcome(n,m);
	for (i = 0; i < n; i++)
		for (j = 0; j < m;j++)
		{
			outcome.M[i][j] = M[i][j] + M1.M[i][j];
		}
	return outcome;
}

Matrix Matrix::operator-(Matrix& M1)
{
	if (m != M1.m || n != M1.n)
	{
		cout << "两矩阵维数不一致，无法进行相减" << endl;
		exit(1);
	}
	int i = 0, j = 0;
	Matrix outcome(n, m);
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			outcome.M[i][j] = M[i][j] - M1.M[i][j];
		}
	return outcome;
}

Matrix Matrix::operator/(double d)
{
	if (d == 0)
	{
		cout << "除数为0，出错退出"<<endl;
		exit(0);
	}
	else
	{
		Matrix TEMP(n, m);
		int i, j;
		for (i = 0; i < n; i++)
			for (j = 0; j< m; j++)
				TEMP.M[i][j] = M[i][j] / d;
		return TEMP;
	}
}

Matrix 	Matrix::Givens(double x, double y)
{
	Matrix TEMP(2,3);
	double c,r,s,t,i,j;
	if (y==0) 
	{
		c = 1;
		s = 0;
		r = x;
	}
	else
	{
		if (abs(y) > abs(x))
		{
			t = -x / y;
			s = sqrt(1 + t*t);
			r = -y*s;
			s = 1 / s;
			c = s*t;
		}
		else
		{
			t = -y / x;
			c = sqrt(1 + t*t);
			r = x*c;
			c = 1 / c;
			s = c*t;
		}

	}
	TEMP.M[0][0] = c;
	TEMP.M[0][1] = s;
	TEMP.M[0][2] = r;
	TEMP,M[1][1] = c;
	TEMP.M[1][0] = -s;
	TEMP.M[1][2] = 0;
	return TEMP;


}

Matrix 	Matrix::operator*(Matrix& M1)

{
	if (m != M1.n )
	{
		cout << "前矩阵列数不等于后矩阵行数，无法进行相乘" << endl;
		exit(1);
	}
	int i = 0, j = 0,k=0;
	Matrix outcome(n, M1.m);

	for (i = 0; i < n; i++)
		for (j = 0; j < M1.m; j++)
		{
			for (k = 0; k < m; k++)
				outcome.M[i][j] += M[i][k] * M1.M[k][j];
		}

		return outcome;
}

Matrix Matrix::operator*(double d)
{
		Matrix TEMP(n, m);
		int i, j;
		for (i = 0; i < n; i++)
			for (j = 0; j < m; j++)
				TEMP.M[i][j] = M[i][j] * d;
		return TEMP;
	
}



void Matrix::operator=(Matrix& M1)
{
	int i, j;
	m = M1.m;
	n = M1.n;
	for (i = 0; i < n;i++)
	{
		for (j = 0; j < m; j++)
		{
			M[i][j] = M1.M[i][j];
		}
	}
//	return M1;

}

bool Matrix::operator==(Matrix& M1)
{
	if (n != M1.n || m != M1.m) return 0;
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			if (M[i][j] != M1.M[i][j]) return 0;
		}
	return 1;
}





void Matrix::DET()
{
	if (m != n)
	{
		cout << "对象不是方阵，无法求其行列式" << endl;
	}
	else
	{
		Matrix TEMP(n,n);
		TEMP = *this;
		int i, j, k, l,s=1;
		double	temp = 0.0;
		for (i = 0; i < n; i++)
		{
			if (TEMP.M[i][i] != 0)
				for (j = i + 1; j < n; j++)
				{
					temp = TEMP.M[j][i] / TEMP.M[i][i];
					for (k = 0; k < n; k++)
					{
						TEMP.M[j][k] -= temp*TEMP.M[i][k];
					}
				}
			else
			{
		/*		cout <<"-----------"<< endl;
				show();
				cout << "-----------" << endl;
		*/
				for (k = i + 1; k <= n; k++)
				{
					if (k == n)
					{
						det = 0;
						 system("pause");
					}
					else  if (TEMP.M[k][i] != 0)
					{
						for (l = 0; l < n; l++)
						{
							temp = TEMP.M[i][l];
							TEMP.M[i][l] = TEMP.M[k][l];
							TEMP.M[k][l] = temp;
						}
						s = s*(-1);
						i--;
						break;
					}			
				}
			}
		}
		det = s;
		for (i = 0; i < n; i++)	det = det*TEMP.M[i][i];
	}
}


bool Matrix::CouldLU()
{
	if (n != m)
	{
		cout << "该矩阵不是方阵，无法进行LU分解" << endl;
		return 0;

	}
	else
	{
		int i,j,k;
		for (i = 1; i < n - 1; i++)
		{
			Matrix MM(i, i);
			for (j = 0; j < i; j++)
				for (k = 0; k < i;k++)
				{
						MM.M[j][k] = M[j][k];
			    }
				
			MM.DET();
			if (MM.det == 0) 
			{
				cout << "该方阵第"<<i<<"阶子行列式为0，无法进行LU分解" << endl;
				return 0;
			}

		}
	}
	return 1;

}

Matrix Matrix::Trans()
{
	Matrix TEMP(m,n);
	int i, j;
	for (i = 0; i < n;i++)
		for (j = 0; j < m; j++)
		{
			TEMP.M[j][i] = M[i][j];
		}
	return TEMP;
}


Matrix* Matrix::LU()
{
	if (!CouldLU())
	{
		cout << "该矩阵无法被LU分解" << endl;
		exit(1);
	}
	else   
	{
		LARGE_INTEGER m_nFreq;
		LARGE_INTEGER m_nBeginTime;
		LARGE_INTEGER nEndTime;
		QueryPerformanceFrequency(&m_nFreq); // 获取时钟周期  
		QueryPerformanceCounter(&m_nBeginTime); // 获取时钟计数  
		Matrix TEMP(n, n);
		TEMP = *this;
		Matrix L(n, n), U(n, n);
		//Matrix LU[2] = { L, U };

		int i, j, k;
		for (k = 0; k < n - 1; k++)
		{
			if (fabs(TEMP.M[k][k] + 1.0) == 1.0)
			{
				cout << "\n分解失败！"<<endl;
				exit(1);
			}
			for (i = k + 1; i < n; i++) TEMP.M[i][k] = TEMP.M[i][k] / TEMP.M[k][k];
			for (i = k + 1; i < n ;i++)
			{
				for (j = k + 1; j < n ; j++)
				{
					TEMP.M[i][j] = TEMP.M[i][j] - TEMP.M[i][k] * TEMP.M[k][j];
				}			
			}
		}

		for (i = 0; i < n ; i++)
		{
			for (j = 0; j < i; j++)
			{
				L.M[i][j] = TEMP.M[i][j];
				U.M[i][j] = 0.0;
			}
			L.M[i][i] = 1.0;
			U.M[i][i] = TEMP.M[i][i];
			for (j = i + 1; j < n ; j++)
			{
				L.M[i][j] = 0.0;
				U.M[i][j] = TEMP.M[i][j];
			}

		}
		QueryPerformanceCounter(&nEndTime);
		cout <<"LU时间花销为："<< (double)(nEndTime.QuadPart - m_nBeginTime.QuadPart) * 1000 / m_nFreq.QuadPart << endl;
	/*	cout << "--------------LLLLLLLLLLL--------------------" << endl;

		L.show();
		cout << "--------------UUUUUUUUUUUUU--------------------" << endl;
		U.show();

		cout << "----------------------------------" << endl;
		*/
		Matrix LU[2] = { L, U };
		return LU;

	}

}

bool Matrix::IsSymmetry()
{
	if (n != m) return 0;
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			if (M[i][j] != M[j][i])
				return 0;
		}
	return 1;
}


bool Matrix::IsPD()
{
	if (n != m)  return 0;
	int i, j, k;
	for (i = 1; i < n; i++)
	{
		Matrix MM(i, i);
		for (j = 0; j < i; j++)
			for (k = 0; k < i; k++)
			{
				MM.M[j][k] = M[j][k];
			}
/*		cout << "========矩阵==========" << endl;
		MM.show();
		cout << "--------行列式------------" << endl;
*/
		MM.DET();
		cout << MM.det << endl;
		if (MM.det <= 0)
		{
			cout << "该方阵第" << i << "阶子行列式不为正，不是正定矩阵" << endl;
			return 0;
		}

	}
	return 1;
}



bool Matrix::CouldQR()
{
	if (n < m)
	{
		cout << "QR分解失败" << endl;
		return 0;
	}
	else if (n == m)
	{
		DET();
		if (det == 0) return 0;
		else return 1;
	}
	return 1;
}

Matrix* Matrix::QR()
{
	if (!CouldQR())
	{
		return 0;
	}
	else
	{
	//	CouldQR();
	//	DET();
	//	cout << det << endl;
		LARGE_INTEGER m_nFreq;
		LARGE_INTEGER m_nBeginTime;
		LARGE_INTEGER nEndTime;
		QueryPerformanceFrequency(&m_nFreq); // 获取时钟周期  
		QueryPerformanceCounter(&m_nBeginTime); // 获取时钟计数  

		Matrix* TEMP=new Matrix(n, m);
		TEMP =this;
		Matrix Q(n,n), R(n,m);
		int i, j, k, mm, jj;
		double u,  alpha, w, t;
		for (i = 0; i <n; i++)
			for (j = 0; j < n; j++)
			{
				Q.M[i][j] = 0.0;
				if (i == j) Q.M[i][j] = 1.0;
			}
		mm = m;
		if (n == m) mm = n - 1;
		for (k = 0; k < mm; k++)
		{
			u = 0.0;
			for (i = k; i < n; i++)
			{
				w = fabs(TEMP->M[i][k]);
				if (w > u) u = w;
			}
			alpha = 0.0;
			for (i = k; i <n; i++)
			{
				t = TEMP->M[i][k] / u;
				alpha = alpha + t*t;
			}
			if (TEMP->M[k][k] > 0.0) u = -u;
			alpha = u*sqrt(alpha);
			if (fabs(alpha) + 1.0 == 1.0)
			{
				cout << "\n		QR分解失败" << endl;
				exit(1);
			}
			u = sqrt(2.0*alpha*(alpha - TEMP->M[k][k]));
			if (u+1.0!=1.0)
			{
				TEMP->M[k][k] = (TEMP->M[k][k] - alpha) / u;
				for (i = k + 1; i < n; i++) TEMP->M[i][k] = TEMP->M[i][k] / u;
				for (j = 0; j < n; j++)
				{
					t = 0.0;
					for (jj = k; jj < n; jj++)
						t = t + TEMP->M[jj][k] * Q.M[jj][j];
					for (i = k; i < n; i++)
						Q.M[i][j] = Q.M[i][j] - 2.0*t*TEMP->M[i][k];
				}
				for (j = k + 1; j < m; j++)
				{
					t = 0;
					for (jj = k; jj < n; jj++)
						t = t + TEMP->M[jj][k] * TEMP->M[jj][j];
					for (i = k; i < n; i++)
						TEMP->M[i][j] = TEMP->M[i][j] - 2.0*t*TEMP->M[i][k];
				}
				TEMP->M[k][k] = alpha;
				for (i = k + 1; i < n; i++) TEMP->M[i][k] = 0.0;
			}
		}
		for (i = 0; i < n - 1; i++)
			for (j = i + 1; j < n; j++)
			{
				t = Q.M[i][j];
				Q.M[i][j] = Q.M[j][i];
				Q.M[j][i] = t;
			}
		QueryPerformanceCounter(&nEndTime);

		cout << "===========================" << endl;
		cout << "QR的时间花销为："<<(double)(nEndTime.QuadPart - m_nBeginTime.QuadPart) * 1000 / m_nFreq.QuadPart << endl;
		cout << "===========================" << endl;

		Q.show();
		cout << "---------------------------" << endl;
		TEMP->show();
		cout << "+++++++++++++++++++++++++++" << endl;
		Matrix QR[2] = {Q, *TEMP};
		return QR;
	}
}


bool Matrix::CouldCholesky()
{
	if (!IsSymmetry()) return 0;
	if (!IsPD())      return 0;
	return 1;
}





Matrix Matrix::Cholesky()
{
	if (!CouldCholesky()||M[0][0]<=0)
	{
		cout << "该矩阵无法被Cholesky分解" << endl;
		exit(1) ;
	}

	LARGE_INTEGER m_nFreq;
	LARGE_INTEGER m_nBeginTime;
	LARGE_INTEGER nEndTime;
	QueryPerformanceFrequency(&m_nFreq); // 获取时钟周期  

	QueryPerformanceCounter(&m_nBeginTime); // 获取时钟计数  
	Matrix L(n, n),LT(n,n);
	L = *this;
	//L.show();

	int i,j,k;
	L.M[0][0] = sqrt(L.M[0][0]);
	for (i = 1; i < n; i++) L.M[i][0] = L.M[i][0] / L.M[0][0];
	for (j = 1; j < n; j++)
	{
		for (k = 0; k < j; k++)
		{
			L.M[j][j] = L.M[j][j]- L.M[j][k] * L.M[j][k];
		}
		if (L.M[j][j] <= 0)
		{
			cout << "Cholesky分解失败" << endl;
			exit(1);
		}
		L.M[j][j] = sqrt(L.M[j][j]);
		cout << L.M[j][j] << endl;
		for (i = j + 1; i < n; i++)
		{
			for (k = 0; k < j; k++)
			{
				L.M[i][j] = L.M[i][j] - L.M[i][k] * L.M[j][k];
			}
			L.M[i][j] = L.M[i][j] / L.M[j][j];
		}
	}
	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++)
			L.M[i][j] = 0.0;
	QueryPerformanceCounter(&nEndTime);
	cout <<"Cholesky的时间花销为："<< (double)(nEndTime.QuadPart - m_nBeginTime.QuadPart) * 1000 / m_nFreq.QuadPart << endl;
    LT = L.Trans();
/*	cout << "-----------LLLLLLLL--------------" << endl;
	L.show();
	cout << "-----------LTLTLTLT--------------" << endl;
	LT.show();
	*/
	return L;
}


/*
Matrix* Matrix::FG()
{
	Matrix TEMP = *this;
	int i, j,k;
	double temp = 0;
	for (i = 0; i < n; i++)
	{
		if (TEMP.M[i][i] != 0)
		{
			for (j = i+1; j < n; j++)
			{
				temp = TEMP.M[j][i] / TEMP.M[i][i];
				for (k = 0; k < m; k++)
					TEMP.M[j][k] = TEMP.M[j][k] - TEMP.M[i][k] * temp;					
			}
			}
		for ()
	}
}
*/

Matrix Matrix::GetLine(int m)
{

	if (m>=this->m) 
	{
		cout << "列序号超出该矩阵索引范围，操作失败" << endl;
		exit(0);
	}
	else {
		Matrix TEMP(n, 1);

		int i;
		for (i = 0; i < n; i++)
			TEMP.M[i][0] = M[i][m];
		return TEMP;
	}

}

Matrix Matrix::GetLine(int a, int b)  //取矩阵第a到第b列
{
	if (b<m && a>=0 && b>=a )
	{
		Matrix TEMP(n, b - a + 1);
		int i, j;
		for (i = 0; i < n; i++)
			for (j = 0; j < b - a + 1; j++)
				TEMP.M[i][j] = M[i][a + j];
		return TEMP;
	}
	else
	{
		cout << "索引超出范围" << endl;
		exit(0);
	}

}




Matrix Matrix::GetRow(int N) //取矩阵第N行
{
	if (N >= n)
	{
		cout << "列序号超出该矩阵索引范围，操作失败" << endl;
		exit(0);
	}
	else {
		Matrix TEMP(1, m);

		int i;
		for (i = 0; i < m; i++)
			TEMP.M[0][i] = M[N][i];
		return TEMP;
	}

}

Matrix Matrix::GetRow(int a, int b)  //取矩阵第a到第b行
{
	if (b<n && a >= 0 && b >= a)
	{
		Matrix TEMP(b-a+1,m);
		int i, j;
		for (i = 0; i < b-a+1; i++)
			for (j = 0; j < m; j++)
				TEMP.M[i][j] = M[a+i][ j];
		return TEMP;
	}
	else
	{
		cout << "索引超出范围" << endl;
		exit(0);
	}

}




double Matrix::GetNorm()   //获得范数平方
{
	double temp = 0;
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			temp += M[i][j] * M[i][j];
	temp = sqrt(temp);
	return temp;
}

Matrix Matrix::Householder()
{
	if (m != 1)
	{
		cout << "该矩阵不是列向量，无法进行household变换" << endl;
		exit(0);
	}
	Matrix  Z(n,1) ;
	Z.M[0][0] = 1;
	Matrix  U(n,1);
	//U = TEMP-Z.ShuCheng(this->GetNorm());
	U = *this - Z*(this->GetNorm());
	U = U / (U.GetNorm());
	Matrix I(n, n),H(n,n);
	int i, j;
	for (i = 0; i < n; i++) I.M[i][i] = 1;
	H = I - U*U.Trans() * 2;
	return H;
}



Matrix Matrix::TwoDiag()
{
	if (n < m)
	{
		exit(0);
	}
	else
	{
		Matrix B(m, m);
		int nn, mm,k;
		Matrix v = this->GetLine(0);
		Matrix A = this->GetLine(1, m - 1);
		for (k=0; k < m-2; k++)
		{
			B.M[k][k] = v.GetNorm();
			Matrix* P = new Matrix(v.n, v.n);
			*P = v.Householder();
			cout << "--------Pv-------------" << endl;
			(*P*v).show();
			cout << "----------P----------" << endl;
			P->show();
			cout << "----------A----------" << endl;
			A.show();
			Matrix PA = (*P)*A;
			cout <<"----------B----------"<< endl;
			B.M[k][k + 1] = PA.GetRow(0).GetNorm();
			B.show();
			cout << "---------PA---------" << endl;
			PA.show();
			Matrix AH = PA.GetRow(1, PA.n-1) * PA.GetRow(0).Trans().Householder();  //
			cout << "----------AH-------------" << endl;
			AH.show();
			cout << "----------v--------------" << endl;
			v = AH.GetLine(0);
			v.show();
			cout << "-------------------------" << endl;
			if (AH.m>1) 	A = AH.GetLine(1, AH.m - 1);
		}
		Matrix P = v.Householder();
		cout << "----------P----------" << endl;
		P.show();
		cout << "----------Pv----------" << endl;
		(P*v).show();
		B.M[k][k] = v.GetNorm();
		cout << "----------A----------" << endl;
		A.show();
		Matrix PA = P*A;
		cout << "---------PA---------" << endl;
		PA.show();
		v = PA.GetRow(1,PA.n-1);
		B.M[k][k + 1] = PA.M[0][0];
		cout << "----------B----------" << endl;
		B.show();
		k++;
		cout << "k=" << k << endl;
		cout << "----------v--------------" << endl;
		v.show();
		B.M[k ][k] = v.GetNorm();	
		int i;
		for (i = 0; i < B.n; i++)
		{
			cout << B.M[i][i] << endl;
			if (B.M[i][i] < 0.0000001)
			{
				B = B.GetLine(0, i - 1).GetRow(0, i - 1);
				cout << "----------B----------" << endl;
				return B;
			}
			B.show();
    	}
		return B;
	}	  
}

Matrix*  Matrix::On_Side_Jacobi()
{
	int iterations = 0;
	Matrix V(m,m),G=*this;
	int i, j,k;
	for (i = 0; i < m; i++)   V.M[i][i] = 1;   
	i = 0;
	j = 0;
	double t = 0, cs = 0, sn = 0, l = 0, tmp = 0;
	bool flag = 1;
	while (iterations < 1000)
	{
		iterations++;
		flag = 1;
		double a = 0, b = 0, c = 0;
		for (i = 0; i < m - 1; i++)
		{
			
			for (j = i + 1; j < m; j++)
			{
				Matrix GG = (this->Trans())*(*this);
				
				for (k = 0; k < n; k++)
				{
					a = a + GG.M[k][i] * GG.M[k][i];
					b = b + GG.M[k][j] * GG.M[k][j];
					c=c+GG.M[k][i] * GG.M[k][j];
				}

				if (abs(c) / sqrt(a*b) > 0.001)  flag = 0;

				l = (b - a) / (2 * c);
				t = sign(l) / (abs(l) + sqrt(1 + l*l));
				cs = 1 / sqrt(1 + t*t);
				sn = cs*t;

				for (k = 0; k <n; k++)
				{
					tmp = M[k][i];
					G.M[k][i] = cs*tmp - sn*G.M[k][j];
					G.M[k][j] = sn*tmp + cs*G.M[k][j];
				}
				for (k = 0; k < m; k++)
				{
					tmp = V.M[k][i];
					V.M[k][i] = cs*tmp - sn*V.M[k][j];
					V.M[k][j] = sn*tmp + cs*V.M[k][j];
				}
			}
		}
		if (flag)			break;
	}
	Matrix E(n, m);
	for (i = 0; i < m; i++)
	{
		E.M[i][i] = G.GetLine(i).GetNorm();		
	}
	G.normalize();
	cout << "------------------------iterations-------------------------"<<endl;
	cout << iterations << endl;
	cout << "-----------------------EEEEEEEEEEEEEEEE------------------------" << endl;

	E.show();
	cout << "-----------------------VVVVVVVVVVVVVV------------------------" << endl;

	V.show();

	cout << "-----------------------GGGGGGGGGGGGG------------------------" << endl;

	G.show();


	Matrix Final[3] = {E,G,V};
	return Final;	
}
