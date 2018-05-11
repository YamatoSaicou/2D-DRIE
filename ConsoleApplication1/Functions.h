#pragma once
#include<iostream>
#include<math.h>
#define nHEIGHT 50
#define nWIDTH 100
#define Mask 2
#define Air 3
#define Surface 1
#define N 1000
#define k0  1.68*pow(10, -15)
#define Fr  6 * pow(10, 14)
#define Ee  0.1 * pow(10, -19)
#define kb  1.38065*pow(10, -23)
#define Ts  300
#define c  2 * pow(10, -23)
#define Yast  0.0032*pow(10,8)
#define Ji  4 * pow(10,16)
#define kD 0.8
#define shigema 0.2
//取k0=1.68*10（-15），Fr=5*10^14,Ee=10^-20,Kb=1.38*10^23,Ts=300. 离子辅助刻蚀率YAST Ji离子流量
using namespace std;
int nMASKTHICK;
int nWINLEFT, nWINRIGHT;
cCell aCA[nHEIGHT][nWIDTH];  //nHEIGHT、nWIDTH分别为元胞系统的高和宽
double f(double seta,double y,double x) 
{
	return (cos(fabs(atan(y / x))-seta))*exp(-(seta-1.57)*(seta-1.57)/(2* shigema*shigema))/(shigema*sqrt(2*3.14));
}


void initialSystem()
{

	cin >> nMASKTHICK;
	cout << "请输入窗口的左起始位置\n";
	cin >> nWINLEFT;
	cout << "请输入窗口的右起始位置\n";
	cin >> nWINRIGHT;
	
	for (int y = 0; y <= nMASKTHICK-1; y++)  // nMASKTHICK为掩膜厚度
	{
		for (int x = 0; x < nWIDTH; x++)
		{
			aCA[y][x].nSpecie = Mask;  aCA[y][x].dAmount = 0;
		}
	}  //将掩膜厚度以内的元胞设置为Mask
	for (int y = 0; y <= nMASKTHICK; y++)
	{
		for (int x = nWINLEFT + 1; x < nWINRIGHT; x++) // 窗口参数
		{
			aCA[y][x].nSpecie = Air;   aCA[y][x].dAmount = 0;  //将窗口位置的元胞设置为Air
		}
	}
	for (int x = nWINLEFT + 1; x < nWINRIGHT; x++) // 窗口参数
	{
		aCA[nMASKTHICK + 1][x].nSpecie = Surface;   aCA[nMASKTHICK + 1][x].dAmount = Full; //将窗口位置下面一层的元胞设置为Surface
	}
	cout << "元胞系统初始化完成，刻蚀开始\n";
}
void XUNHUAN()
{

	for (int y = nMASKTHICK; y < nHEIGHT; y++)  // 求出表面元胞法向量（全部的细胞扫描一次的方法来找到表面元胞细胞）
	{
		for (int x = 0; x < nWIDTH; x++)
		{
			if (aCA[y][x].nSpecie == Surface)
			{
				double p[1][2] = { 0,0 }; //P1 2 3 4 5 6 7 8是表面元胞相邻八个元胞的表面法向量，a 1-8是对应的气体贡献度
				double p1[1][2] = { 0,0 },a1=0;
				double p2[1][2] = { 0,0 },a2=0;
				double p3[1][2] = { 0,0 },a3=0;
				double p4[1][2] = { 0,0 },a4=0;
				double p5[1][2] = { 0,0 }, a5 = 0;
				double p6[1][2] = { 0,0 }, a6 = 0;
				double p7[1][2] = { 0,0 }, a7 = 0;
				double p8[1][2] = { 0,0 }, a8 = 0; 
				if ((aCA[y + 1][x].nSpecie) == Air) { p1[0][1] = { 1 }, a1 = 1; }
				if ((aCA[y-1][x].nSpecie) == Air) { p2[0][1] = { -1 }, a2 = 1; }
				if ((aCA[y][x-1].nSpecie) == Air) { p3[0][0] = { -1 }, a3 = 1; }
				if ((aCA[y][x+1].nSpecie) == Air) { p4[0][0] = { 1 }, a4 = 1; }
				if ((aCA[y + 1][x + 1].nSpecie) == Air) { p5[0][1] = { 0.707 },p5[0][0]= 0.707, a5 = 0.707; }
				if ((aCA[y - 1][x - 1].nSpecie) == Air) { p6[0][1] = { -0.707 },p6[0][0]=-0.707, a6 = 0.707; }
				if ((aCA[y - 1][x + 1].nSpecie) == Air) { p7[0][1] = { -0.707 }, p7[0][0] = { 0.707 }, a7 = 0.707; }
				if ((aCA[y + 1][x - 1].nSpecie) == Air) { p8[0][0] = { 0.707 }, p8[0][1] = { -0.707 }, a8 = 0.707; }
				p[0][0] = { p1[0][0] + p2[0][0] + p3[0][0] + p4[0][0] + p5[0][0] + p6[0][0] + p7[0][0] + p8[0][0] };
				p[0][1] = { p1[0][1] + p2[0][1] + p3[0][1] + p4[0][1]+ p5[0][1] + p6[0][1] + p7[0][1] + p8[0][1] };//得到法向量的两个分量
			/*	cout << "坐标是是(" << x << "," << y << ")" << "的表面元胞法向量是" << "(" << p[0][0] << "," << p[0][1] << ")" << "\n"; */
				//接下来计算刻蚀速率
				double b1, b2;
				 b1 = (double)p[0][0];
				 b2 = (double)p[0][1];
				double v1 = 0.0000; double v2 = 0;
				double seta1 = 0, seta2 = 0;
				double m=0;
				m= k0 * Fr*exp(-Ee / (kb*Ts));
				v1 = m*(a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8);   //化学刻蚀速率 
				if (x < nWINLEFT)
				{
					seta1 = fabs(atan(((double)y - nMASKTHICK) / ((double)x - nWINLEFT)));
					seta2 = fabs(atan(((double)y - 0) / ((double)x - nWINRIGHT)));
				}
				else if (x == nWINLEFT)
				{
					seta1 = 3.14 / 2;
					seta2 = fabs(atan(((double)y - 0) / ((double)x - nWINRIGHT)));
				}
				else if (x == nWINRIGHT)
				{
					seta1 = fabs(atan(((double)y - 0) / ((double)x - nWINLEFT)));
					seta2 = 3.14 / 2;
				}
				else if (x>nWINRIGHT)
				{
					seta1 = fabs(atan(((double)y) / ((double)x - nWINLEFT)));
					seta2 = fabs(atan(((double)y - nMASKTHICK) / ((double)x - nWINRIGHT)));
				}
				else if (nWINLEFT< x <nWINRIGHT)
				{
					seta1 = 3.14-fabs(atan(((double)y - 0) / ((double)x - nWINLEFT)));
					seta2 = fabs(atan(((double)y - 0) / ((double)x - nWINRIGHT)));
				}
				double sum = 0.0;                               //用科斯特方法求积分
				double seta = fabs((seta1 - seta2));
				double gaps = fabs((seta1 - seta2)) / double(1000);  //每个间隔的长度
				double h = gaps / 2.0;
				for (int i = 0; i < 1000; i++)
				{
						sum += (h / 45.0) * fabs(((7.0*f(i * gaps, b2, b1))) +
							fabs(32.0*f(i * gaps + 0.25*gaps, b2, b1)) +
							fabs(12.0*f( i * gaps + 0.5*gaps, b2, b1)) +
							fabs(32.0*f( i * gaps + 0.75*gaps, b2, b1)) +
							fabs(7.0*f((i + 1)*gaps, b2, b1)));
				}
			/*	cout << "积分结果是" << sum << " " << atan(b2 / b1)<<" "<<b1<<" "<<b2; */
				double Ci=0;
				Ci = c * Yast *Ji;
				v2 = Ci*sum;
				double v = v1+ v2;
			/*	cout << "坐标是" << x << "," << y << "seta1和seta2的差是" << (fabs(seta2 - seta1) / 3.14) * 180 << "度" << seta1 / 3.14 * 180 << "度" << seta2 / 3.14 * 180 << "v1是"<<v1<<"v2是"<<v2<<"\n";
			/*	cout << "坐标是是(" << x << "," << y << ")" << "的元胞刻蚀速率是" << v << "\n"; */
				aCA[y][x].dAmount = aCA[y][x].dAmount - v;
			/*	cout << "坐标是是(" << x << "," << y << ")" << "的元胞剩余含量是" << aCA[y][x].dAmount << "\n"; */
				if (aCA[y][x].dAmount <= 0)
				{
					aCA[y][x].nSpecie = Air;
					aCA[y][x].dAmount = 0;
					/* cout << "坐标是是(" << x << "," << y << ")" << "的元胞被置为空气元胞\n"; */
					if (aCA[y + 1][x].nSpecie == Si)
					{
						aCA[y + 1][x].nSpecie = Surface;
						/*cout << "将其上方的元胞置为表面元胞" << "\n";*/
					}
					if (aCA[y - 1][x].nSpecie == Si)
					{
						aCA[y - 1][x].nSpecie = Surface;
						/*cout << "将其下方的元胞置为表面元胞" << "\n";*/
					}
					if (aCA[y][x - 1].nSpecie == Si)
					{
						aCA[y][x - 1].nSpecie = Surface;
						/*cout << "将其左方的元胞置为表面元胞" << "\n";*/
					}
					if (aCA[y][x + 1].nSpecie == Si)
					{
						aCA[y][x + 1].nSpecie = Surface;
						/*cout << "将其右方的元胞置为表面元胞" << "\n";*/
					}
				}
				if (aCA[y][x].dAmount > 0)
				{
					aCA[y][x].dAmount = aCA[y][x].dAmount + v1 *kD;
					/*cout << "坐标是是(" << x << "," << y << ")" << "的沉积后含量是" << aCA[y][x].dAmount << "\n";*/
				}

			}
		}
	} 




}