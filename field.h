#pragma once

class field
{
public:
	virtual ~field(void);
	 field(void);
	virtual void calc_field()=0;
	//virtual double get_field(double x1, double x2);
};
