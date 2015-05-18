
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <bitset>

inline void printBitVector( std::vector<bool> bitvec )
{
	for ( unsigned int i = 0; i<bitvec.size(); ++i)
	{
		//std::cout << i << " = " << (bitvec[i] ? '1' : '0') << " : ";
		std::cout << (bitvec[i] ? '1' : '0');
	}	
	std::cout << std::endl;
}

template<unsigned int T>
inline void printBitVector( std::bitset<T> bitvec )
{
	for ( unsigned int i = 0; i<bitvec.size(); ++i)
	{
		//std::cout << i << " = " << (bitvec[i] ? '1' : '0') << " : ";
		std::cout << (bitvec[i] ? '1' : '0');
	}	
	std::cout << std::endl;
}

template<typename T>
inline void printVector( std::vector<T> vec )
{
	for ( unsigned int i = 0; i<vec.size(); ++i)
	{
		std::cout << i << " = " << vec[i] << " : ";		
	}	
	std::cout << std::endl;
}

template<typename T>
inline void printVector( std::vector<T> vec, unsigned int n )
{
	for ( unsigned int i = 0; i<n; ++i)
	{
		std::cout << i << " = " << vec[i] << " : ";		
	}	
	std::cout << std::endl;
}

template<typename T>
inline void printVector( T* vec, unsigned int n )
{
	for ( unsigned int i = 0; i<n; ++i)
	{
		std::cout << i << " = " << vec[i] << " : ";		
	}	
	std::cout << std::endl;
}

template<typename T>
inline void printVector( std::vector<T> vec, unsigned int rows, unsigned int cols )
{
	for ( unsigned int i = 0; i<rows; ++i)
	{
		for ( unsigned int j = 0; j<cols; ++j)
		{
			std::cout << vec[j + i*cols] << " : ";		
		}
		std::cout << std::endl;
	}	
}

template<typename T>
inline void printVector( T* vec, unsigned int rows, unsigned int cols )
{
	for ( unsigned int i = 0; i<rows; ++i)
	{
		for ( unsigned int j = 0; j<cols; ++j)
		{
			std::cout << vec[j + i*cols] << " : ";
		}
		std::cout << std::endl;
	}	
}

template<typename T>
inline void writeFile( std::string filename, T* vec, unsigned int rows, unsigned int cols )
{
	std::ofstream myfile;
	myfile.open ( filename );

	for ( unsigned int i = 0; i<rows; ++i)
	{
		for ( unsigned int j = 0; j<cols; ++j)
		{
			 myfile << std::setprecision(8) << std::setw(15) << vec[j + i*cols] << " : ";			
		}
		  myfile << "\n";
	}	

	myfile.close();
}

template<typename T>
inline void writeFile( std::string filename, std::vector<T> vec, unsigned int rows, unsigned int cols )
{
	std::ofstream myfile;
	myfile.open ( filename );

	for ( unsigned int i = 0; i<rows; ++i)
	{
		for ( unsigned int j = 0; j<cols; ++j)
		{
			 myfile << std::setprecision(8) << std::setw(15) << vec[j + i*cols] << " : ";
		}
		  myfile << "\n";
	}	

	myfile.close();
}

template<typename T>
inline void writeFileTabbed( std::string filename, std::vector<T> vec, unsigned int rows, unsigned int cols )
{
	std::ofstream myfile;
	myfile.open ( filename );

	for ( unsigned int i = 0; i<rows; ++i)
	{
		for ( unsigned int j = 0; j<cols; ++j)
		{
			 myfile << vec[j + i*cols] << '\t';
		}
		  myfile << "\n";
	}	

	myfile.close();
}

template<typename T>
inline void writeFileTabbed( std::string filename, T* vec, unsigned int rows, unsigned int cols )
{
	std::ofstream myfile;
	myfile.open ( filename );

	for ( unsigned int i = 0; i<rows; ++i)
	{
		for ( unsigned int j = 0; j<cols; ++j)
		{
			 myfile << vec[j + i*cols] << '\t';
		}
		  myfile << "\n";
	}	

	myfile.close();
}

template<typename T>
inline void writeFileTabbed( std::string filename, std::vector<T> vec )
{
	std::ofstream myfile;
	myfile.open ( filename );

	for ( unsigned int i = 0; i<vec.size(); ++i)
	{
		myfile << vec[i] << '\n';		
	}

	myfile.close();
}

template<typename T>
inline void writeFileTabbed( std::string filename, std::vector<T> vec, unsigned int n )
{
	std::ofstream myfile;
	myfile.open ( filename );

	for ( unsigned int i = 0; i<n; ++i)
	{
		myfile << vec[i] << '\n';		
	}

	myfile.close();
}

template<typename T>
inline void readFile( std::string filename, std::vector<T>& vec )
{
	std::ifstream myfile( filename );

	std::string line;
	int i = 0;
	while (getline(myfile, line))
	{
		std::stringstream  lineStream(line);
		std::string  cell;
		while (getline(lineStream, cell, ','))
		{
			double ds = atof(cell.c_str());
			vec[i] = ds;
			i++;
		}
	}

	myfile.close();
}

template<typename T>
inline bool compareVectors( std::vector<T>& vec1, std::vector<T>& vec2 )
{
	bool flag = true;
	for ( int i = 0; i < vec1.size(); i++ )
	{	
		if ( abs( vec1[i] - vec2[i] ) > 0.001 ) flag = false;
	}
	return flag;
}