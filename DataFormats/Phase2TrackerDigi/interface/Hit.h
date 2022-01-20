#ifndef HIT_H
#define HIT_H

class Hit {
	private:
		int row;
		int col;
		int adc;

	public:
		Hit(int row_num, int col_num, int adc_num);

		int row() const {
			return row;
		};

		int col() const {
			return col;
		};

		int adc() const {
			return adc;
		};
};

#endif
