
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//Тип состояния. Плоско-напряжённое (длина >> ширины), Плоско-деформированное.
enum stressType {stress, strain};

void ATBAMult(const std::vector<std::vector <double> > B, const std::vector<std::vector <double> >& D, std::vector<std::vector <double> >& result)
{
	for(size_t i =0; i<result.size(); i++)
		for(size_t j=0; j<result[i].size(); j++)
			result[i][j]=0;

	for (size_t i = 0; i<result.size(); i++)
	{
		for (size_t l = 0; l < result.size(); l++)
		{
			double temp2=0;
			for (size_t j = 0; j < D.size(); j++)
			{
				double temp=0;
				for (size_t k = 0; k < D.size(); k++)
				{
					temp+=D[j][k]*B[k][l];
				}
				temp2+=B[j][i]*temp;
			}
			result[i][l] = temp2;
		}
	}
}

double Ddot(std::vector<double>& l, std::vector<double>& r)
{
	double result=0;
	if (l.size()!=r.size())
		return result;
	for(size_t i=0; i<l.size();i++)
		result+=l[i]*r[i];
	return result;

}

void MatrVectMul(std::vector<double>& q, const std::vector<std::vector<double>>& A, const std::vector<double> p)
{
	for(size_t i=0; i<A.size(); i++)
	{
		q[i]=0;
		for(size_t j=0; j<A.size(); j++)
		{
			q[i]+=A[i][j]*p[j];
		}
	}
}

void Daxpy(std::vector<double>& result, double beta, std::vector<double>& second, double alpha)
{
	for(size_t i = 0; i<result.size(); i++)
		result[i]=beta*result[i]+alpha*second[i];
}


class Solver
{
	std::vector<std::vector<double> > A;
	std::vector<double> b;
	std::vector<double> x;
	size_t max_iter;
	double tol;
	Solver(){max_iter=0; tol = 0.0;}
public:
	Solver(std::vector<std::vector<double>>& TangentMatrix, std::vector<double>& resist_Force, size_t max_iter = 1000, double tol = 1e-6)
	: A(TangentMatrix), b(resist_Force), max_iter(max_iter), tol(tol){}

	std::vector<double> CG()
	{
		std::cout<<"Begin Solver CG on CPU from CSR matrix"<<std::endl;

		std::cout<<"A:"<<std::endl;
		for(auto i: A)
		{
			for(auto j:i)
				std::cout<<j<<" ";
			std::cout<<std::endl;
		}

		std::cout<<"b:"<<std::endl;
		for(auto i: b)
		{
			std::cout<<i<<" ";
		}
		std::cout<<std::endl;

		x.resize(b.size(), 0);
		std::vector<double> r(b), p(b), q(b.size());
		double alpha=0, beta=0, nrmr, ro, nrmb=Ddot(b,b);
		ro=nrmr=nrmb;
		nrmr=nrmb=sqrt(nrmb);

		std::cout<<nrmb<<std::endl;

		size_t k=0;
		while ((nrmr/nrmb > tol)&&(k<max_iter))
		{
			MatrVectMul(q,A,p);
			alpha = Ddot(q,p);

			alpha = ro/alpha;
			Daxpy(x, 1, p, alpha);
			Daxpy(r, 1, q, -alpha);
			beta = ro;
			ro = Ddot(r,r);
			nrmr=sqrt(ro);
			beta = ro/beta;
			Daxpy(p, beta, r, 1);
			k++;
			if (!(k%100)) std::cout<<k<<": nrmr = "<<nrmr<<std::endl;
		}

		std::cout<<"End Solver CG on CPU from CSR matrix \n"
				 <<"Drop tolerance: "<<sqrt(nrmr/nrmb)<<"\n"
				 <<"Number of iterations: "<<k<<"\n"
				 <<"Solve norma: "<<sqrt(Ddot(x,x))<<std::endl;
		return x;
	}
};



class Dot
{
	double x,y;
	size_t global_number;
public:
	Dot (double x=0, double y=0, size_t global_number=0): x(x),y(y), global_number(global_number){};
	double get_X() {return x;}
	double get_Y() {return y;}
	size_t get_Number() {return global_number;}
	double operator[] (int i)
	{
		if (i==0)
			return x;
		if (i==1)
			return y;
		return 0;
	}
	Dot operator- (Dot& second)
	{
		return Dot (x - second.x, y - second.y);
	}
	void Set (double x, double y, size_t global_number) {this->x=x; this->y=y; this->global_number = global_number; }
};

class Material
{
	stressType _stressType;
	double E; //модуль Юнга
	double nu; //коэффициент Пуассона
	std::vector<std::vector<double> > D;

public:
	Material(double E, double nu, stressType _stressType = stress): E(E), nu(nu), _stressType(_stressType)
	{
		D.resize(3);
		D[0].resize(3); D[1].resize(3); D[2].resize(3);
		double d1, d2, b;
		if (_stressType == stress)
		{
			// plain stress плоское напряжение - пластины (продольное сечение стержня - длина>>ширины)
			d1 = E / (1.0 - nu * nu);
			d2 =  E * 0.5/(1.0 + nu);
			b  = d1 * nu;

		}
		else if (_stressType == strain)
		{
			// plain strain - плоская деформация, стержни
			d1 = E * (1-nu) / ((1.0 + nu)*(1-2*nu));
			d2 = E * 0.5 / (1 + nu);
			b = E * nu / ((1.0 + nu) * (1.0 - 2*nu));
		}
		D[0][0] = d1; D[0][1] =  b; D[0][2] =  0;
		D[1][0] =  b; D[1][1] = d1; D[1][2] =  0;
		D[2][0] =  0; D[2][1] =  0; D[2][2] = d2;
	}
	std::vector<std::vector<double> > & getMaterialMatrix()
	{
		return D;
	}
};

class Quad
{
	std::vector<Dot> Node;
	Material* _Material;
	std::vector<std::vector<double> >TangentMatrix;
	std::vector<double> resist_Force;
	std::vector<double> mass_Matrix;
	size_t _size;
	size_t number_of_nodes;
	double load_force;
public:
	Quad(): _size(0), number_of_nodes(0), _Material(NULL), load_force(0){};
	~Quad()
	{}
	size_t get_Number(size_t i)
	{
		return Node[i].get_Number();
	}
	size_t get_Number_of_Nodes()
	{
		return number_of_nodes;
	}
	void FormHex(std::vector<Dot> Node, Material* _Material)
	{

		this->Node = Node;
		this->_Material = _Material;
		_size = 8;
		number_of_nodes = 4;
		TangentMatrix.resize(_size, std::vector<double>(_size, 0));
		resist_Force.resize(_size);
		mass_Matrix.resize(3*3);
		load_force=0;
	}
	void set_Load(double load)
	{
		load_force = load;
	}
	std::vector<std::vector<double> > FormStiffnesMatrix()
	{
		const short OrderOfQuadrature = 2;
		// узлы интегрирования для квадратур Гаусса
		double qGauss[OrderOfQuadrature];
		// веса для квадратур Гаусса
		double Weights[OrderOfQuadrature];
		qGauss[0] = 0.577350; qGauss[1] = -0.577350;
		Weights[0] = Weights[1] = 1;

		std::vector<std::vector <double> > BtDB (_size, std::vector<double> (_size,0));

		for(size_t i=0; i<TangentMatrix.size(); i++)
			for(size_t j=0; j<TangentMatrix[i].size(); j++)
				TangentMatrix[i][j]=0;

		for (unsigned i = 0; i < OrderOfQuadrature; i++)
			for (unsigned int j = 0; j < OrderOfQuadrature; j++)
			{
				// узлы интегрирования Гаусса
				double r = qGauss[i];
				double s = qGauss[j];

				double detJacoby = 0;
				// вычисляем матрицу производных B от функции форм
				std::vector<std::vector<double> > B, D;
				B = ComputeBn(r, s, detJacoby);
				//матрица материала
				D=_Material->getMaterialMatrix();

				ATBAMult(B,D,BtDB);

				double V = (fabs(detJacoby) * Weights[j] * Weights[i]);

				for(size_t i=0; i<BtDB.size(); i++)
					for(size_t j=0; j<BtDB[i].size(); j++)
						BtDB[i][j] *= V;

				for(size_t i=0; i<TangentMatrix.size(); i++)
					for(size_t j=0; j<TangentMatrix[i].size(); j++)
						TangentMatrix[i][j]+=BtDB[i][j];
			}
		return TangentMatrix;
	}

	// сформировать локальный вектор правых частей
	std::vector<double> ResistingForce()
	{
		int i, j, k;
		int numEdges = 4; // 4 ребра
		for(size_t i=0; i<resist_Force.size(); i++)
			resist_Force[i]=0;

		if (fabs(load_force)>1e-6)
		{
			double half_length = sqrt((Node[3].get_X()-Node[2].get_X())*(Node[3].get_X()-Node[2].get_X())+(Node[3].get_Y()-Node[2].get_Y())*(Node[3].get_Y()-Node[2].get_Y()))*0.5;
			resist_Force[5] += load_force*half_length;
			resist_Force[7] += load_force*half_length;
		}

		return resist_Force;
	}

	double ComputeInverseJacoby (const double r, const double s, std::vector<Dot>& deltaCoor, std::vector<std::vector <double> >& InverseJacoby)
	{

		std::vector <std::vector<double> > Jacoby(2, std::vector<double> (2,0));
		for (int i = 0; i < 2; i++)
		{
			Jacoby[0][i] = 0.25 * ((1 - s) * deltaCoor[0][i] + (1 + s) * deltaCoor[1][i]);
			Jacoby[1][i] = 0.25 * ((1 - r) * deltaCoor[2][i] + (1 + r) * deltaCoor[3][i]);
		}

		double det = Jacoby[0][0] * Jacoby[1][1] - Jacoby[0][1] * Jacoby[1][0];

		InverseJacoby[0][0] = Jacoby[1][1] / det;
		InverseJacoby[1][1] = Jacoby[0][0] / det;
		InverseJacoby[0][1] = -Jacoby[0][1] / det;
		InverseJacoby[1][0] = -Jacoby[1][0] / det;

		return det;
	}

	std::vector<std::vector<double> > ComputeBn (const double r, const double s, double& detJacoby)
	{
		std::vector<std::vector <double> > InverseJacoby (2, std::vector<double> (2,0));
		std::vector<Dot> deltaCoor(4);
		deltaCoor[0] = Node[1] - Node[0];
		deltaCoor[1] = Node[2] - Node[3];
		deltaCoor[2] = Node[3] - Node[0];
		deltaCoor[3] = Node[2] - Node[1];

		detJacoby = ComputeInverseJacoby(r, s, deltaCoor, InverseJacoby);
		// производные по глобальным координатам x,y
		std::vector<double> dNx(4);
		std::vector<double> dNy(4);

		// производные по локальным координатам r,s
		std::vector<double> dNr(4);
		std::vector<double> dNs(4);

		//вычисляем производные от билинейных функций форм по локальным координатам r,s
		Compute_dNrs(r,s,dNr,dNs);

		//вычисляем производные от функций форм по глобальным координатам x,y
		for (int k = 0; k < 4; k++)
		{
			dNx[k] = InverseJacoby[0][0] * (dNr)[k] + InverseJacoby[0][1] * (dNs)[k];
			dNy[k] = InverseJacoby[1][0] * (dNr)[k] + InverseJacoby[1][1] * (dNs)[k];
		}


		std::vector<std::vector<double> > B (3, std::vector<double> (8,0));

		//заполняем матрицу производных от функции форм B
		ComputeB(dNx, dNy, B);

		return B;
	}

	//вычисление производных от функций форм по локальным координатам r,s
	void Compute_dNrs(const double r, const double s, std::vector<double>& dNr, std::vector<double>& dNs) const
	{
		(dNr)[0] = 0.25 * (s - 1);	(dNs)[0] = 0.25 * (r - 1);
		(dNr)[1] = 0.25 * (1 - s);	(dNs)[1] = 0.25 * (-1 - r);
		(dNr)[2] = 0.25 * (s + 1);	(dNs)[2] = 0.25 * (1 + r);
		(dNr)[3] = 0.25 * (-s - 1);	(dNs)[3] = 0.25 * (1 - r);
	}

	// формирование матрицы B производных от функции форм
	void ComputeB(const std::vector<double>& dNx, const std::vector<double> dNy, std::vector<std::vector<double> > &B) const
	{
		for (int i = 0; i < 4; i++)
		{
			B[0][2 * i] = B[2][2 * i + 1] = dNx[i];
			B[1][2 * i + 1] = B[2][ 2 * i] = dNy[i];
		}
	}

};

class mesh
{
	std::vector<Quad> Element;
	std::vector<double> coordinates;
	size_t Number_of_points;
	std::vector<std::vector<double> > TangentMatrix;
	std::vector<double> resist_Force;
	std::vector<std::pair<unsigned,double>> ConstraintNodes;
	Material* _material;
	std::vector<double> displacement;
public:
	mesh(): _material(NULL){};

	std::vector<std::vector<double> > FormStiffnesMatrix ()
	{

		for(size_t i=0; i<TangentMatrix.size(); i++)
			for(size_t j=0; j<TangentMatrix[i].size(); j++)
			{
				TangentMatrix[i][j]=0;
			}

		for(size_t i=0; i<Element.size(); i++)
		{
			std::vector<std::vector<double> > LocalStiffnesMatrix = Element[i].FormStiffnesMatrix();

			std::cout<<"Local matrix "<<i<<":"<<std::endl;
			for(auto k: LocalStiffnesMatrix)
			{
				for(auto j:k)
					std::cout<<j<<" ";
				std::cout<<std::endl;
			}

			for(size_t j=0; j<Element[i].get_Number_of_Nodes(); j++)
			{
				for(size_t k=0; k<Element[i].get_Number_of_Nodes(); k++)
				{
					TangentMatrix[2*Element[i].get_Number(j)][2*Element[i].get_Number(k)] += LocalStiffnesMatrix[2*j][2*k];
					TangentMatrix[2*Element[i].get_Number(j)+1][2*Element[i].get_Number(k)+1] += LocalStiffnesMatrix[2*j+1][2*k+1];
				}
			}
		}

		return TangentMatrix;
	}

	std::vector<double> ResistingForce()
	{
		for(size_t i=0; i<resist_Force.size(); i++)
			resist_Force[i]=0;
		for(size_t i=0; i<Element.size(); i++)
		{
			std::vector<double> LocalResistForce = Element[i].ResistingForce();

			for(size_t j=0; j<Element[i].get_Number_of_Nodes(); j++)
			{
				resist_Force[2*Element[i].get_Number(j)] += LocalResistForce[2*j];
				resist_Force[2*Element[i].get_Number(j)+1] += LocalResistForce[2*j+1];
			}
		}

		return resist_Force;
	}

	void Constraint()
	{
		for(size_t i=0; i<ConstraintNodes.size(); i++)
		{
			unsigned number = ConstraintNodes[i].first;
			double t = ConstraintNodes[i].second;
			TangentMatrix[2*number][2*number] = 1.0;
			resist_Force[2*number] = t;
			for(size_t j=0; j<TangentMatrix.size(); j++)
			{
				if(j!=2*number)
				{
					TangentMatrix[2*number][j]=0;
					resist_Force[j]-=TangentMatrix[j][2*number]*t;
					TangentMatrix[j][2*number]=0;
				}
			}
			TangentMatrix[2*number+1][2*number+1] = 1;
			resist_Force[2*number+1] = t;
			for(size_t j=0; j<TangentMatrix.size(); j++)
			{
				if(j!=2*number+1)
				{
					TangentMatrix[2*number+1][j]=0;
					resist_Force[j]-=TangentMatrix[j][2*number+1]*t;
					TangentMatrix[j][2*number+1]=0;
				}
			}
		}
	}

	void SetMesh(std::vector<double>& coordinate, unsigned nX, unsigned nY, std::vector<std::pair<unsigned,double> > constraint, std::vector<unsigned>& load, double load_force,
				 double E, double Nu, stressType _stressType = stress)
	{

		_material = new Material(E, Nu, _stressType);
		Element.resize(0);
		ConstraintNodes = constraint;
		Number_of_points = coordinate.size()/2;
		coordinates = coordinate;

		for(unsigned i = 0; i<nY-1; i++)
		{
			for(unsigned j = 0; j<nX-1; j++)
			{
				std::vector<Dot> _quad;
				//unsigned number_of_nodes = 4;
				_quad.push_back(Dot(coordinate[2*(nX*i + j)],coordinate[2*(nX*i + j) + 1], nX*i + j));
				_quad.push_back(Dot(coordinate[2*(nX*i + j + 1)],coordinate[2*(nX*i + j + 1) + 1], nX*i + j + 1));
				_quad.push_back(Dot(coordinate[2*(nX*(i+1) + j)],coordinate[2*(nX*(i+1) + j) + 1], nX*(i+1) + j));
				_quad.push_back(Dot(coordinate[2*(nX*(i+1) + j + 1)],coordinate[2*(nX*(i+1) + j + 1) + 1], nX*(i+1) + j + 1));
				Quad _el;
				_el.FormHex(_quad, _material);
				Element.push_back(_el);
			}
		}

		for(unsigned i=0; i<load.size(); i++)
			Element[load[i]].set_Load(load_force);

		TangentMatrix.resize(coordinate.size(), std::vector<double>(coordinate.size(), 0));
		resist_Force.resize(coordinate.size(), 0);
	}

	void Solve()
	{
		std::cout<<"FormStiffnesMatrix"<<std::endl;
		FormStiffnesMatrix();
		std::cout<<"ResistingForce"<<std::endl;
		ResistingForce();

		std::cout<<"Constraint"<<std::endl;
		Constraint();

		Solver S(TangentMatrix, resist_Force);
		displacement = S.CG();
	}

	void Print()
	{
		std::ofstream out("out.vtu");
		out<<"<?xml version=\"1.0\"?>"<<std::endl
		   <<"<VTKFile type=\"UnstructuredGrid\" vesion=\"0.1\" byte_order=\"LittleEndian\">"<<std::endl
		   <<"  <UnstructuredGrid>"<<std::endl
		   <<"    <Piece NumberOfPoints=\""<<Number_of_points<<"\" NumberOfCells=\""<<Element.size()<<"\">"<<std::endl
		   <<"      <PointData Scalars=\"scalars\">"<<std::endl
		   <<"        <DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" Format=\"ascii\">"<<std::endl;
		for(unsigned i=0; i<Number_of_points; i++)
			out<<displacement[2*i]<<" "<<displacement[2*i+1]<<" "<<0<<" ";
		out<<std::endl;
		out<<"        </DataArray>"<<std::endl
		   <<"      </PointData>"<<std::endl
		   <<"      <Points>"<<std::endl
		   <<"        <DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" Format=\"ascii\">"<<std::endl;
		for(unsigned i=0; i<Number_of_points; i++)
			out<<coordinates[2*i]<<" "<<coordinates[2*i+1]<<" "<<0<<" ";
		out<<std::endl;
		out<<"        </DataArray>"<<std::endl
		   <<"      </Points>"<<std::endl
		   <<"      <Cells>"<<std::endl
		   <<"        <DataArray type=\"Int64\" Name=\"connectivity\" Format=\"ascii\">"<<std::endl;
		for(auto i:Element)
			out<<i.get_Number(0)<<" "<<i.get_Number(1)<<" "<<i.get_Number(2)<<" "<<i.get_Number(3)<<" ";
		out<<std::endl;
		out<<"        </DataArray>"<<std::endl
		   <<"        <DataArray type=\"Int64\" Name=\"offsets\" Format=\"ascii\">"<<std::endl;
		unsigned k=4;
		for(unsigned i=0; i<Element.size(); i++)
		{
			out<<k<<" ";
			k+=4;
		}
		out<<std::endl;
		out<<"        </DataArray>"<<std::endl
		   <<"        <DataArray type=\"Int64\" Name=\"types\" Format=\"ascii\">"<<std::endl;
		for(unsigned i=0; i<Element.size(); i++)
		{
			out<<9<<" ";
		}
		out<<std::endl;
		out<<"        </DataArray>"<<std::endl
		   <<"      </Cells>"<<std::endl
		   <<"    </Piece>"<<std::endl
		   <<"</UnstructuredGrid>"<<std::endl
		   <<"</VTKFile>"<<std::endl;
		out.close();
	}
	void PrintDisplacement()
	{
		std::ofstream out("displacement.txt");
		out<<"const point = ["<<std::endl;
		for(unsigned i=0; i<Number_of_points; i++)
			out<<"{ x: "<<displacement[2*i]<<", y:"<<displacement[2*i+1]<<"},"<<std::endl;
		out<<"];"<<std::endl;
		out.close();
	}
};


int main()
{
	std::cout<<"Begin programm"<<std::endl;
	std::vector<double> coor;
//	double X = 10.0, Y = 3.0;
//	size_t Nx = 4;
//	size_t Ny = 3;
	double X = 12.0, Y=5.0;
	size_t Nx = 20, Ny = 10;
	size_t NumberOfNode = (Nx+1)*(Ny+1);
	const size_t NoD=2;

	double dx = X/Nx;
	double dy = Y/Ny;
	double x=0.0, y=0.0;
	std::cout<<"Form Mesh"<<std::endl;
	while(y<=Y)
	{
		while (x<=X)
		{
			coor.push_back(x);
			coor.push_back(y);
			x+=dx;
		}
		y+=dy;
		x=0.0;
	}
	
	std::vector<std::pair<unsigned,double>> constraint = {{20,0.0} ,{41,0.0},{62,0.0},{83,0.0},{104,0.0},{125,0.0},{146,0.0},{167,0.0},{188,0.0},{209,0.0},{230,0.0}};



	std::vector<unsigned> load = {180,181, 182, 183, 184, 185, 186, 187, 188, 189};

	double load_force = 10;
	double Nu=0.34;
	double E = 70;
	stressType _stressType = stress;

	std::cout<<coor.size()<<std::endl;

	mesh Test;
	std::cout<<"Form Element Mesh"<<std::endl;
	Test.SetMesh(coor, Nx+1, Ny+1, constraint, load, load_force, E, Nu, _stressType);
	std::cout<<"Solving.."<<std::endl;
	Test.Solve();
	std::cout<<"Print"<<std::endl;
	Test.Print();
	Test.PrintDisplacement();

	std::cout<<"End programm"<<std::endl;

}

