#ifndef _GRID_H_
#define _GRID_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "global.h"
#include "bnd_cond.h"

const unsigned int CELL_FLAG_GOOD	= 0x000000;
const unsigned int CELL_FLAG_BAD	= 0x000001;
const unsigned int CELL_FLAG_LIM	= 0x000002;



class Cell 
{
public:
	Cell(): nodesInd(NULL), edgesInd(NULL), nCount(0), eCount(0) {};
	~Cell();

	int   nCount;
	int   eCount;
	int*  nodesInd;
	int*  edgesInd;
	int	  neigh[3];
	int	  type;
	char  typeName[64];
	double  S;
	Point c; // ����� ������
	double HX;
	double HY;
	unsigned int flag;

	friend class Grid;
};

class Edge 
{
public:
	Edge(): c(NULL), cCount(0) {};
	~Edge();

	int      n1;        // ���� � ������
	int      n2;        // ���� � �����
	int      c1;        // ������ �����
	int      c2;        // ������ ������, ������� �� �1 � �2
	Vector   n;         // ������� � �����
	double     l;         // ����� ����� 
	double     cnl1;       // ����� ���� �������� �������� �� ������ ����� �� ������ ������� ������
	double     cnl2;       // ����� ���� �������� �������� �� ������ ����� �� ������ ������� ������
	int      cCount;    // ���������� ����� �� �����
	Point*   c;         // ����� �� �����
	int      type;      // ��� ����� (�����., ������.)
	char	typeName[64];
	CFDBoundary * bnd = NULL;
	friend class Grid;
public:
	static const int TYPE_INNER		= 0;  
	static const int TYPE_INLET		= 1;  
	static const int TYPE_OUTLET	= 2;  
	static const int TYPE_WALL		= 3;  

	static const int TYPE_NAMED = 0x100;
};

class Grid 
{
public:
	Grid(): nodes(NULL), cells(NULL), edges(NULL),
		nCount(0), cCount(0), eCount(0) {};
	~Grid();

	void initFromFiles(char* fName);
	inline Point& getNode(int i) { return nodes[i]; };
	inline Cell&  getCell(int i) { return cells[i]; };
	inline Edge&  getEdge(int i) { return edges[i]; };
	void replaceEdges(int if1, int if2);
	virtual void reorderEdges();
	int findEdge(int n1, int n2);

	Point* nodes;
	Cell*  cells;
	Edge*  edges;
	int    nCount;
	int    cCount;
	int    eCount;
};

inline double _max_(double a, double b)
{
	if (a>b) { return a; } else {return b; }
}

inline double _min_(double a, double b)
{
	if (a<b) { return a; } else {return b; }
}

inline double _max_(double a, double b, double c)
{
	return _max_(a, _max_(b, c));
}

inline double _min_(double a, double b, double c)
{
	return _min_(a, _min_(b, c));
}

inline double _sqr_(double a)
{
	return a*a;
}

#endif