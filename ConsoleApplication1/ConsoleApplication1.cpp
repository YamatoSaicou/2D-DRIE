#include "stdafx.h"
#include<iostream>
#include <fstream>
#include<cmath>
#include"Classes.h"
#include"Functions.h"

using namespace std;
int main()
{
	cout<<"请输入掩膜厚度\n";
	initialSystem();
	for (int i = 0; i < 100; i++)
	{
		XUNHUAN();
		cout << "第"<<i+1<<"次刻蚀/沉积循环完成\n";
	}
	ofstream outfile;
	string b = "result.txt";
	outfile.open(b);
	for (int y = 0; y < nHEIGHT; y++)
	{
		for (int x = 0; x < nWIDTH; x++)
		{
			outfile << aCA[y][x].nSpecie;
			outfile << " ";
		}
		outfile << "\n";
	}
    return 0;
}

