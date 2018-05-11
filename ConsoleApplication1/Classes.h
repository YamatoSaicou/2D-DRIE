#pragma once
#define Si 4
#define Mask 2
#define Air 3
#define Full 1
class cCell {
public:
	double dAmount;  //元胞的种类参数和容量参数
	int nSpecie;
	cCell()
	{
		nSpecie = Si; dAmount = Full;
	}  //默认种类为Si，默认容量为Full
	cCell(int a, double b)
	{
		nSpecie = a; dAmount = b;
	}  //构造函数
};
