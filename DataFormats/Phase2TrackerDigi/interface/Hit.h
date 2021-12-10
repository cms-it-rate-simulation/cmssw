#ifndef HIT_H
#define HIT_H

class Hit {
	public:
		int row;
		int col;
		int adc;

		Hit(int row_num, int col_num, int adc_num);
};

#endif
